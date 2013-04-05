(ns bcbio.variation.report
  "Parse and provide detailed information from GATKReport outputs."
  (:use [ordered.map :only [ordered-map]]
        [clojure.math.combinatorics :only [cartesian-product]]
        [bcbio.variation.variantcontext :only [parse-vcf get-vcf-retriever
                                               variants-in-region
                                               get-vcf-iterator]]
        [bcbio.variation.callable :only [get-callable-checker is-callable?]]
        [bcbio.variation.evaluate :only [organize-gatk-report-table]]
        [bcbio.variation.metrics :only [ml-on-vcf-metrics passes-filter? nonref-passes-filter?]])
  (:require [clojure.string :as string]
            [clojure.data.csv :as csv]
            [doric.core :as doric]
            [lonocloud.synthread :as ->]
            [bcbio.variation.grade :as grade]))

(defn concordance-report-metrics
  "Retrieve high level concordance metrics from GATK VariantEval report."
  [sample in-file]
  (letfn [(sample-in-row? [x]
            (and (= (:Sample x) sample)
                 (= (:Novelty x) "all")
                 (= (:Filter x) "called")))]
    (->> (organize-gatk-report-table in-file "GenotypeConcordance" sample-in-row?)
         (map (fn [x] [(keyword (:variable x)) (:value x)]))
         (into {}))))

(defn count-variants
  "Count variants that pass an optional checker function."
  [f ref-file check?]
  (with-open [vcf-iter (get-vcf-iterator f ref-file)]
    (count (filter check? (parse-vcf vcf-iter)))))

(defn shared-discordant
  "Provide metrics to distinguish types of discordance in a comparison.
  These identify variants which differ due to being missing in one variant
  call versus calls present in both with different genotypes. It also pulls
  out variants in diploids that differ due to het/hom calls."
  [file1 file2 ref-file]
  (with-open [file1-iter (get-vcf-iterator file1 ref-file)
              vcf-retriever (get-vcf-retriever ref-file file2)]
    (reduce (fn [coll vc]
              (let [other-vcs (variants-in-region vcf-retriever
                                                  (:chr vc) (:start vc) (:end vc))
                    vc-type (if-not (empty? other-vcs) :total :unique)]
                (-> coll
                    (assoc vc-type (inc (get coll vc-type)))
                    (assoc :hethom ((if (and (= :total vc-type)
                                             (grade/hethom-discordant? vc other-vcs))
                                      inc identity)
                                    (get coll :hethom))))))
            {:total 0 :unique 0 :hethom 0}
            (parse-vcf file1-iter))))

(defn nocoverage-count
  "Calculate count of variant in input file without coverage in the comparison."
  [in-vcf ref-file compare-kw compared]
  (let [out-dir (get-in compared [:dir :prep] (get-in compared [:dir :out]))
        align-file (get-in compared [compare-kw :align]
                           (get-in compared [:exp :align]))]
    (when-not (nil? align-file)
      (with-open [call-source (get-callable-checker align-file (-> compared :exp :ref)
                                                    :out-dir out-dir)]
        (count-variants in-vcf ref-file
                        #(and (passes-filter? %)
                              (not (is-callable? call-source (:chr %) (:start %) (:end %)))))))))

(defn get-summary-level
  "Retrieve expected summary level from configuration"
  [config]
  (letfn [(level-from-string [x]
            (case (when-not (nil? x) (string/lower-case x))
              "quick" :quick
              "full" :full
              :standard))
          (get-string-level [config to-try]
            (loop [cur-try to-try]
              (if (empty? cur-try) nil
                  (let [cur-level (get-in config (first cur-try))]
                    (if-not (nil? cur-level) cur-level
                            (recur (rest cur-try)))))))]
    (let [to-check (cartesian-product [:exp :c1 :c2 :call] [:summary-level])]
      (level-from-string (get-string-level config to-check)))))

(defn get-sv-metrics
  "Retrieve structural variation metrics from SV concordance files."
  [compared ref-file]
  (when-let [finfo (seq (->> (:c-files compared)
                             (drop-while #(not= (first %) :sv-concordant))
                             (take 3)))]
    (reduce (fn [coll [kw vcf-file]]
              (assoc coll kw
                     (ordered-map
                      :total (count-variants vcf-file ref-file passes-filter?))))
            (ordered-map) finfo)))

(defn- all-vrn-counts
  "Provide counts in a file, split by type of variation."
  [fname cmp-kw compared]
  (letfn [(vrn-type-passes-filter? [vrn-type]
            (fn [vc]
              (and (passes-filter? vc)
                   (contains? vrn-type (:type vc)))))]
    (let [sum-level (get-summary-level compared)
          ref-file (get-in compared [:exp :ref])
          base {:total (count-variants fname ref-file passes-filter?)}]
      (if (= sum-level :quick) base
          (assoc base
            :nocoverage (nocoverage-count fname ref-file cmp-kw compared)
            :snp (count-variants fname ref-file
                                 (vrn-type-passes-filter? #{"SNP"}))
            :indel (count-variants fname ref-file
                                   (vrn-type-passes-filter? #{"INDEL"})))))))

(defn- grading-discordance-metrics
  "Prepare grading based metrics to break down reasons for discordant calls"
  [vcf-file]
  (println vcf-file))

(defn- add-discordance-metrics
  "Add detailed discordance metrics to the summary information."
  [m compared]
  (let [ref-file (get-in compared [:exp :ref])
        c-files (-> compared :c-files vals)]
    (-> m
        (assoc :discordant1 (all-vrn-counts (second c-files) :c2 compared))
        (->/when (> (count c-files) 2)
          (assoc :discordant2 (all-vrn-counts (nth c-files 2) :c1 compared))
          (assoc :discordant_both (apply shared-discordant (conj (vec (take 2 (rest c-files)))
                                                                 ref-file)))))))

(defn- add-concordance-metrics
  "Add detailed concordance metrics to summary information."
  [m compared]
  (let [ref-file (get-in compared [:exp :ref])
        c-files (-> compared :c-files vals)]
    (-> m
        (assoc :genotype_concordance (-> compared :metrics :percent_overall_genotype_concordance))
        (assoc :callable_concordance (-> compared :callable-metrics
                                         :percent_overall_genotype_concordance))
        (assoc :nonref_discrepency (-> compared :metrics :percent_non_reference_discrepancy_rate))
        (assoc :nonref_sensitivity (-> compared :metrics :percent_non_reference_sensitivity))
        (assoc :concordant (all-vrn-counts (first c-files) nil compared))
        (assoc :nonref_concordant (count-variants (first c-files) ref-file
                                                  nonref-passes-filter?)))))

(defn top-level-metrics
  "Provide one-line summary of similarity metrics for a VCF comparison."
  [compared]
  (let [sum-level (get-summary-level compared)
        ref-file (get-in compared [:exp :ref])
        c-files (-> compared :c-files vals)]
    (-> (ordered-map)
        (assoc :sample (-> compared :exp :sample))
        (assoc :ftypes (take 3 (-> compared :c-files keys)))
        (add-concordance-metrics compared)
        (add-discordance-metrics compared)
        (->/when-let [sv-metrics (get-sv-metrics compared ref-file)]
          (assoc :sv sv-metrics))
        (->/when (= sum-level :full)
          (assoc :ml_metrics (ml-on-vcf-metrics ref-file (take 2 c-files)))))))

(defn calc-accuracy
  "Calculate an overall accuracy score from input metrics.
  The accuracy logic is:
  (#correctly aligned bases - (1*(simple substitutions and indels) +
                               2*(larger errors))
   / #correctly aligned bases)"
  [metrics error-items]
  (letfn [(get-penalty [[error-type call-type]]
            (case call-type
              :snp 1
              :indel 2
              :sv 2))]
    (let [error-items (cartesian-product error-items [:snp :indel :sv])
          error-score (apply + (map #(* (get-in metrics % 0) (get-penalty %)) error-items))
          total-bases (get-in metrics [:total-bases :compared] 1)]
      (float
       (* 100.0 (/ (- total-bases error-score) total-bases))))))

(defn prep-scoring-table
  "Summary table of high level variables and scoring metrics for comparison."
  [metrics sv-stats]
  (let [to-write (ordered-map :accuracy "Accuracy score"
                              :accuracy-phasing "Accuracy score, including phasing"
                              [:total-bases :percent] "Completeness"
                              [:total-bases :compared] "Total bases scored"
                              [:total-bases :total] "Possible evaluation bases"
                              [:discordant :snp] "Discordant SNPs"
                              [:discordant :indel] "Discordant indels"
                              [:discordant :sv] "Discordant structural variants"
                              [:phasing-error :snp] "Phasing Error SNPs"
                              [:phasing-error :indel] "Phasing Error indels"
                              :haplotype-blocks "Phased haplotype blocks"
                              ;:nonmatch-het-alt "Non-matching heterozygous alternative alleles"
                              )
        sv-metrics (assoc-in metrics [:discordant :sv]
                             (-> sv-stats vals second (get :total 0)))
        s-metrics (-> sv-metrics
                      (assoc :accuracy (calc-accuracy sv-metrics [:discordant]))
                      (assoc :accuracy-phasing (calc-accuracy sv-metrics
                                                              [:discordant :phasing-error])))
        need-percents {:accuracy 6
                       :accuracy-phasing 6
                       [:total-bases :percent] 2}]
    (letfn [(prep-row [[k x]]
              (let [val (if (coll? k) (get-in s-metrics k) (get s-metrics k))]
                {:metric x
                 :value (if (contains? need-percents k)
                          (format (str "%." (get need-percents k) "f") val)
                          val)}))]
      (map prep-row to-write))))

(defn- write-scoring-table
  "Write high level metrics table in readable format."
  [metrics sv-stats wrtr]
  (when-not (or (nil? metrics)
                (nil? (get-in metrics [:total-bases :total])))
    (.write wrtr (str (doric/table [:metric :value] (prep-scoring-table metrics sv-stats))
                      "\n"))))

(defn- write-concordance-metrics
  "Summary table of metrics for assessing the score of a variant comparison."
  [metrics wrtr]
  (letfn [(metrics-info [ftype-i & kvs]
            (if (<= (count (:ftypes metrics)) ftype-i)
              []
              (let [cur-name (name (nth (:ftypes metrics) ftype-i))]
                (apply concat
                       (map (fn [[k v]] [k (str cur-name ": " v)]) (partition 2 kvs))))))]
    (let [to-write (apply ordered-map
                          (concat [:genotype_concordance "Overall genotype concordance"
                                   :callable_concordance "Callable genotype concordance"
                                   :nonref_discrepency "Non-reference discrepancy rate"
                                   :nonref_sensitivity "Non-reference sensitivity"]
                                  (metrics-info 0
                                                [:concordant :total] "total"
                                                :nonref_concordant "non-reference"
                                                [:concordant :snp] "SNPs"
                                                [:concordant :indel] "indels")
                                  (metrics-info 1
                                                [:discordant1 :total] "total"
                                                [:discordant1 :nocoverage] "unique"
                                                [:discordant1 :snp] "SNPs"
                                                [:discordant1 :indel] "indels")
                                  (metrics-info 2
                                                [:discordant2 :total] "total"
                                                [:discordant2 :nocoverage] "unique"
                                                [:discordant2 :snp] "SNPs"
                                                [:discordant2 :indel] "indels")
                                  [[:discordant_both :total] "Shared discordant"
                                   [:discordant_both :hethom] "het/hom discordant"
                                   [:ml_metrics :top-metrics] "Classification metrics"]))]
      (letfn [(get-value [[k metric]]
                (when-let [val (if (coll? k) (get-in metrics k) (get metrics k))]
                  {:metric metric :value val}))]
        (.write wrtr (str (doric/table [:metric :value] (remove nil? (map get-value to-write)))
                          "\n"))))))

(defn- write-sv-metrics
  "Summary table of structural variation comparions."
  [sv-metrics wrtr]
  (letfn [(get-values [[base xs]]
            (map (fn [[inner-kw val]]
                   {:metric (str (name base) ": " (name inner-kw))
                    :value val})
                 xs))]
    (.write wrtr "** Structural variation\n")
    (.write wrtr (str (doric/table [:metric :value]
                                   (->> (map get-values sv-metrics)
                                        flatten
                                        (remove nil?)))
                      "\n"))))

;; ## Classification metrics

(defn- write-classification-metrics
  "Summary table of classification metrics from GATK variant recalibration."
  [cmp-info wrtr]
  (letfn [(get-metric-counts [in-vcf]
            (with-open [vcf-iter (get-vcf-iterator in-vcf (get-in cmp-info [:exp :ref]))]
              (reduce (fn [coll vc]
                        (let [culprit (get-in vc [:attributes "culprit"])]
                          (if (or (nil? culprit) (= (count (:filters vc)) 0)) coll
                              (assoc coll culprit (inc (get coll culprit 0))))))
                      {} (parse-vcf vcf-iter))))
          (get-recal-metrics [in-vcf]
            (sort-by :count >
                     (map (fn [[m c]] {:metric m :count c}) (get-metric-counts in-vcf))))]
    (.write wrtr "** GATK recalibration filter metrics\n")
    (doseq [call (map (partial get cmp-info) [:c1 :c2])]
      (when (= (:mod call) "recal")
        (.write wrtr (str (doric/table [:metric :count]
                                       (get-recal-metrics (:file call)))
                          "\n"))))))

;; ## Top level reports

(defn write-summary-txt
  "Write a summary text file with tables of useful concordance metrics."
  [wtr comparisons]
  (doseq [x comparisons]
    (.write wtr (format "* %s : %s vs %s\n" (get-in x [:exp :sample])
                        (get-in x [:c1 :name]) (get-in x [:c2 :name])))
    (write-scoring-table (:metrics x) (get-in x [:summary :sv]) wtr)
    (write-concordance-metrics (:summary x) wtr)
    (when-let [sv-info (get-in x [:summary :sv])]
      (write-sv-metrics sv-info wtr))
    (when (get-in x [:c1 :mod])
      (write-classification-metrics x wtr))))

(defn- get-summary-csv-vals
  "Retrieve values for CSV output checking for lots of special cases.
    - Display all values for :total lines
    - For :snp and :indel lines, only show values with dictionaries that
      have this detailed info.
    - Always display initial sample and call values (i <= 3)"
  [header cmp cur-type]
  (for [v (map-indexed (fn [i k]
                         (let [v (get cmp k)]
                           (cond
                            (= :type k) (name cur-type)
                            (and (not= :total cur-type)
                                 (> i 3)
                                 (not (map? v))) nil
                            :else v)))
                       header)]
    (if (map? v) (get v cur-type) v)))

(defn write-summary-csv
  "Write a top level summary CSV file with useful concordance metrics."
  [wtr comparisons]
  (doseq [[i [x cmp-orig]] (map-indexed vector (map (juxt identity :summary) comparisons))]
    (let [header (concat [:sample :call1 :call2 :type] (nnext (keys cmp-orig)))
          cmp (-> cmp-orig
                  (dissoc :ftypes)
                  (assoc :call1 (-> x :c1 :name))
                  (assoc :call2 (-> x :c2 :name)))]
      (when (= i 0)
        (csv/write-csv wtr [(map name header)]))
      (csv/write-csv wtr (for [cur-type [:total :snp :indel]]
                           (get-summary-csv-vals header cmp cur-type))))))

(defn write-files-csv
  "Write a CSV file with locations of useful files."
  [wtr comparisons config]
  (csv/write-csv wtr [["call1" "call2" "type" "fname"]])
  (doseq [x comparisons]
    (doseq [[k f] (:c-files x)]
      (csv/write-csv wtr [[(get-in x [:c1 :name]) (get-in x [:c2 :name]) (name k)
                           (string/replace-first f (str (get-in config [:dir :out]) "/") "")]]))))
