(ns bcbio.variation.report
  "Parse and provide detailed information from GATKReport outputs."
  (:import [org.broadinstitute.sting.gatk.report GATKReport])
  (:use [ordered.map :only [ordered-map]]
        [clojure.math.combinatorics :only [cartesian-product]]
        [bcbio.variation.variantcontext :only [parse-vcf get-vcf-retriever
                                               get-vcf-source]]
        [bcbio.variation.callable :only [callable-checker]])
  (:require [doric.core :as doric]))

(defn concordance-report-metrics
  "Retrieve high level concordance metrics from GATK VariantEval report."
  [sample in-file]
  (letfn [(sample-in-row? [x]
            (and (= (:row x) sample)
                 (= (:Sample x) sample)
                 (= (:Novelty x) "all")
                 (= (:Filter x) "called")))]
    (let [table (-> (GATKReport. in-file) (.getTable "GenotypeConcordance.simplifiedStats"))
          cols (rest (.getColumns table))
          headers (map #(-> % (.getColumnName) keyword) cols)]
      (filter sample-in-row?
              (for [i (range (count (.values (first cols))))]
                (zipmap headers
                        (map #(nth (vec (.values %)) i) cols)))))))

(defn- count-variants
  "Count variants that pass an optional checker function."
  [f check?]
  (with-open [vcf-source (get-vcf-source f)]
    (count (filter check? (parse-vcf vcf-source)))))

(defn discordance-metrics
  "Provide metrics to distinguish types of discordance in a comparison.
  These identify variants which differ due to being missing in one variant
  call versus calls present in both with different genotypes."
  [file1 file2]
  (with-open [file2-source (get-vcf-source file2)
              file1-source (get-vcf-source file1)]
    (let [vrn-fetch (get-vcf-retriever file2-source)]
      (reduce (fn [coll vc]
                (let [other-vcs (vrn-fetch (:chr vc) (:start vc) (:end vc))
                      vc-type (if-not (empty? other-vcs) :total :unique)]
                  (assoc coll vc-type (inc (get coll vc-type)))))
              {:total 0 :unique 0}
              (parse-vcf file1-source)))))

(defn nocoverage-count
  "Calculate count of variant in input file without coverage in the comparison."
  [in-vcf compare-kw compared]
  (let [align-file (get-in compared [compare-kw :align]
                           (get-in compared [:exp :align]))]
    (if (nil? align-file)
      ""
      (let [[callable? call-source] (callable-checker align-file (-> compared :exp :ref)
                                                      :out-dir (-> compared :dir :out))
            vc-callable? (fn [vc]
                           (callable? (:chr vc) (:start vc) (:end vc)))]
        (with-open [_ call-source]
          (count-variants in-vcf vc-callable?))))))

(defn top-level-metrics
  "Provide one-line summary of similarity metrics for a VCF comparison."
  [compared]
  (letfn [(passes-filter? [vc]
            (= (count (:filters vc)) 0))
          (nonref-passes-filter? [vc]
            (and (passes-filter? vc)
                 (every? #(contains? #{"HET" "HOM_VAR"} (:type %)) (:genotypes vc))))
          (vrn-type-passes-filter [vrn-type]
            (fn [vc]
              (and (passes-filter? vc)
                   (contains? vrn-type (:type vc)))))
          (all-vrn-counts [fname cmp-kw compared]
            {:total (count-variants fname passes-filter?)
             :nocoverage (nocoverage-count fname cmp-kw compared)
             :snp (count-variants fname (vrn-type-passes-filter #{"SNP"}))
             :indel (count-variants fname (vrn-type-passes-filter #{"INDEL"}))})]
    (ordered-map
     :sample (-> compared :exp :sample)
     :call1 (-> compared :c1 :name)
     :call2 (-> compared :c2 :name)
     :genotype_concordance (-> compared :metrics :percent_overall_genotype_concordance)
     :nonref_discrepency (-> compared :metrics :percent_non_reference_discrepancy_rate)
     :nonref_sensitivity (-> compared :metrics :percent_non_reference_sensitivity)
     :concordant (all-vrn-counts (first (:c-files compared)) nil compared)
     :nonref_concordant (count-variants (first (:c-files compared)) nonref-passes-filter?)
     :discordant1 (all-vrn-counts (second (:c-files compared)) :c2 compared)
     :discordant2 (all-vrn-counts (nth (:c-files compared) 2) :c1 compared)
     :discordant_both (apply discordance-metrics (rest (:c-files compared))))))

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
  (.write wrtr (str (doric/table [:metric :value] (prep-scoring-table metrics))
                    "\n")))

(defn calc-score
  "Calculate scoring for input metric types"
  [type val]
  (if (keyword? type) ""
      (let [sign (if (= :concordant (first type)) 1 -1)]
        (case (second type)
              :snp (-> 1 (* sign) (* val))
              :indel (-> 2 (* sign) (* val))
              0))))

(defn write-concordance-metrics
  "Summary table of metrics for assessing the score of a variant comparison."
  [metrics wrtr]
  (let [to-write (ordered-map :genotype_concordance "Overall genotype concordance"
                              :nonref_discrepency "Non-reference discrepancy rate"
                              :nonref_sensitivity "Non-reference sensitivity"
                              [:concordant :total] "Total concordant"
                              :nonref_concordant "Non-reference concordant count"
                              [:discordant_both :total] "Shared discordant"
                              [:concordant :snp] "Concordant SNPs"
                              [:concordant :indel] "Concordant indels"
                              [:discordant1 :total] (str "Discordant total: " (:call1 metrics))
                              [:discordant1 :nocoverage]  (str "Discordant unique: "
                                                               (:call1 metrics))
                              [:discordant1 :snp] (str "Discordant SNPs: " (:call1 metrics))
                              [:discordant1 :indel] (str "Discordant indels: " (:call1 metrics))
                              [:discordant2 :total] (str "Discordant total: " (:call1 metrics))
                              [:discordant2 :nocoverage]  (str "Discordant unique: "
                                                               (:call2 metrics))
                              [:discordant2 :snp] (str "Discordant SNPs: " (:call2 metrics))
                              [:discordant2 :indel] (str "Discordant indels: " (:call2 metrics)))]
    (letfn [(get-value-and-score [[k metric]]
              (let [val (if (coll? k) (get-in metrics k) (get metrics k))]
                {:metric metric :value val :score (calc-score k val)}))]
      (.write wrtr (str (doric/table [:metric :value]
                                     (map get-value-and-score to-write))
                        "\n")))))
