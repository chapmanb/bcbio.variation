(ns bcbio.variation.report
  "Parse and provide detailed information from GATKReport outputs."
  (:import [org.broadinstitute.sting.gatk.report GATKReport])
  (:use [ordered.map :only [ordered-map]]
        [clojure.math.combinatorics :only [cartesian-product]]
        [bcbio.variation.variantcontext :only [parse-vcf get-vcf-retriever
                                               get-vcf-source]]
        [bcbio.variation.callable :only [callable-checker]]
        [bcbio.variation.metrics :only [ml-on-vcf-metrics passes-filter? nonref-passes-filter?]])
  (:require [doric.core :as doric]
            [clojure.string :as string]))

(defn concordance-report-metrics
  "Retrieve high level concordance metrics from GATK VariantEval report."
  [sample in-file]
  (letfn [(sample-in-row? [x]
            (and (= (:Sample x) sample)
                 (= (:Novelty x) "all")
                 (= (:Filter x) "called")))]
    (let [table (-> (GATKReport. in-file) (.getTable "GenotypeConcordance"))
          cols (rest (.getColumns table))
          headers (map #(-> % (.getColumnName) keyword) cols)]
      (->> (for [i (range (count (.values (first cols))))]
             (zipmap headers
                     (map #(nth (vec (.values %)) i) cols)))
           (filter sample-in-row?)
           (map (fn [x] [(keyword (:variable x)) (:value x)]))
           (into {})))))

(defn- count-variants
  "Count variants that pass an optional checker function."
  [f ref-file check?]
  (with-open [vcf-source (get-vcf-source f ref-file)]
    (count (filter check? (parse-vcf vcf-source)))))

(defn discordance-metrics
  "Provide metrics to distinguish types of discordance in a comparison.
  These identify variants which differ due to being missing in one variant
  call versus calls present in both with different genotypes."
  [file1 file2 ref-file]
  (with-open [file2-source (get-vcf-source file2 ref-file)
              file1-source (get-vcf-source file1 ref-file)]
    (let [vrn-fetch (get-vcf-retriever file2-source)]
      (reduce (fn [coll vc]
                (let [other-vcs (vrn-fetch (:chr vc) (:start vc) (:end vc))
                      vc-type (if-not (empty? other-vcs) :total :unique)]
                  (assoc coll vc-type (inc (get coll vc-type)))))
              {:total 0 :unique 0}
              (parse-vcf file1-source)))))

(defn nocoverage-count
  "Calculate count of variant in input file without coverage in the comparison."
  [in-vcf ref-file compare-kw compared]
  (let [align-file (get-in compared [compare-kw :align]
                           (get-in compared [:exp :align]))]
    (if (nil? align-file)
      ""
      (let [out-dir (get-in compared [:dir :prep] (get-in compared [:dir :out]))
            [callable? call-source] (callable-checker align-file (-> compared :exp :ref)
                                                      :out-dir out-dir)
            vc-callable? (fn [vc]
                           (callable? (:chr vc) (:start vc) (:end vc)))]
        (with-open [_ call-source]
          (count-variants in-vcf ref-file vc-callable?))))))

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
  [finfo ref]
  (reduce (fn [coll [kw vcf-file]]
            (assoc coll kw
                   (ordered-map
                    :total (count-variants vcf-file ref passes-filter?))))
          (ordered-map) finfo))

(defn top-level-metrics
  "Provide one-line summary of similarity metrics for a VCF comparison."
  [compared]
  (let [sum-level (get-summary-level compared)
        ref-file (get-in compared [:exp :ref])]
    (letfn [(vrn-type-passes-filter [vrn-type]
              (fn [vc]
                (and (passes-filter? vc)
                     (contains? vrn-type (:type vc)))))
            (all-vrn-counts [fname cmp-kw compared]
              (let [base {:total (count-variants fname ref-file passes-filter?)}]
                (if (= sum-level :quick) base
                    (assoc base
                      :nocoverage (nocoverage-count fname ref-file cmp-kw compared)
                      :snp (count-variants fname ref-file
                                                (vrn-type-passes-filter #{"SNP"}))
                      :indel (count-variants fname ref-file
                                             (vrn-type-passes-filter #{"INDEL"}))))))]
      (let [c-files (-> compared :c-files vals)
            base
            (ordered-map
             :sample (-> compared :exp :sample)
             :ftypes (take 3 (-> compared :c-files keys))
             :sv (let [xs (->> (:c-files compared)
                               (drop-while #(not= (first %) :sv-concordant))
                               (take 3))]
                   (when-not (empty? xs)
                     (get-sv-metrics xs ref-file)))
             :genotype_concordance (-> compared :metrics :percent_overall_genotype_concordance)
             :callable_concordance (-> compared :callable-metrics
                                       :percent_overall_genotype_concordance)
             :nonref_discrepency (-> compared :metrics :percent_non_reference_discrepancy_rate)
             :nonref_sensitivity (-> compared :metrics :percent_non_reference_sensitivity)
             :concordant (all-vrn-counts (first c-files) nil compared)
             :nonref_concordant (count-variants (first c-files) ref-file
                                                nonref-passes-filter?)
             :discordant1 (all-vrn-counts (second c-files) :c2 compared)
             :discordant2 (when (> (count c-files) 2) (all-vrn-counts (nth c-files 2) :c1 compared))
             :discordant_both (when (> (count c-files) 2)
                                (apply discordance-metrics (conj (vec (take 2 (rest c-files)))
                                                                 ref-file))))]
        (if-not (= sum-level :full) base
            (assoc base
              :ml_metrics (ml-on-vcf-metrics ref-file (take 2 c-files))))))))

(defn calc-accuracy
  "Calculate an overall accuracy score from input metrics.
  The accuracy logic is:
  (#correctly aligned bases / (#correctly aligned bases +
                               1*(simple substitutions and indels) +
                               2*(larger errors)))."
  [metrics]
  (letfn [(get-penalty [[error-type call-type]]
            (case call-type
              :snp 1
              :indel 2))]
    (let [error-items (cartesian-product [:discordant :phasing-error] [:snp :indel])
          error-score (apply + (map #(* (get-in metrics % 0) (get-penalty %)) error-items))
          total-bases (get-in metrics [:total-bases :compared] 1)]
      (float
       (* 100.0 (/ total-bases (+ total-bases error-score)))))))

(defn prep-scoring-table
  "Summary table of high level variables and scoring metrics for comparison."
  [metrics]
  (let [to-write (ordered-map :accuracy "Overall accuracy score"
                              [:total-bases :percent] "Percentage of bases compared"
                              [:total-bases :compared] "Total bases compared"
                              [:total-bases :total] "Possible evaluation bases"
                              [:discordant :snp] "Discordant SNPs"
                              [:discordant :indel] "Discordant indels"
                              [:phasing-error :snp] "Phasing Error SNPs"
                              [:phasing-error :indel] "Phasing Error indels"
                              :haplotype-blocks "Phased haplotype blocks"
                              :nonmatch-het-alt "Non-matching heterozygous alternative alleles")
        s-metrics (assoc metrics :accuracy (calc-accuracy metrics))
        need-percent #{:accuracy [:total-bases :percent]}]
    (letfn [(prep-row [[k x]]
              (let [val (if (coll? k) (get-in s-metrics k) (get s-metrics k))]
                {:metric x
                 :value (if (contains? need-percent k) (format "%.2f" val) val)}))]
      (map prep-row to-write))))

(defn write-scoring-table
  "Write high level metrics table in readable format."
  [metrics wrtr]
  (when-not (or (nil? metrics)
                (nil? (get-in metrics [:total-bases :total])))
    (.write wrtr (str (doric/table [:metric :value] (prep-scoring-table metrics))
                      "\n"))))

(defn write-concordance-metrics
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
                                   [:ml_metrics :top-metrics] "Classification metrics"]))]
      (letfn [(get-value [[k metric]]
                (when-let [val (if (coll? k) (get-in metrics k) (get metrics k))]
                  {:metric metric :value val}))]
        (.write wrtr (str (doric/table [:metric :value] (remove nil? (map get-value to-write)))
                          "\n"))))))

(defn write-sv-metrics
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

(defn write-classification-metrics
  "Summary table of classification metrics from GATK variant recalibration."
  [cmp-info wrtr]
  (letfn [(get-metric-counts [in-vcf]
            (with-open [vcf-source (get-vcf-source in-vcf (get-in cmp-info [:exp :ref]))]
              (reduce (fn [coll vc]
                        (let [culprit (get-in vc [:attributes "culprit"])]
                          (if (or (nil? culprit) (= (count (:filters vc)) 0)) coll
                              (assoc coll culprit (inc (get coll culprit 0))))))
                      {} (parse-vcf vcf-source))))
          (get-recal-metrics [in-vcf]
            (sort-by :count >
                     (map (fn [[m c]] {:metric m :count c}) (get-metric-counts in-vcf))))]
    (.write wrtr "** GATK recalibration filter metrics\n")
    (doseq [call (map (partial get cmp-info) [:c1 :c2])]
      (when (= (:mod call) "recal")
        (.write wrtr (str (doric/table [:metric :count]
                                       (get-recal-metrics (:file call)))
                          "\n"))))))
