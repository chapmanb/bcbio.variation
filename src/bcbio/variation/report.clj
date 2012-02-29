(ns bcbio.variation.report
  "Parse and provide detailed information from GATKReport outputs."
  (:import [org.broadinstitute.sting.gatk.report GATKReport])
  (:use [ordered.map :only [ordered-map]]
        [clojure.math.combinatorics :only [cartesian-product]]
        [bcbio.variation.variantcontext :only [parse-vcf get-vcf-retriever
                                               get-vcf-source]])
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

(defn discordance-metrics
  "Provide metrics to distinguish types of discordance in a comparison.
  These identify variants which differ due to being missing in one variant
  call versus calls present in both with different genotypes."
  [file1 file2]
  (with-open [vcf-source (get-vcf-source file2)]
    (let [vrn-fetch (get-vcf-retriever vcf-source)]
      (reduce (fn [coll vc]
                (let [other-vcs (vrn-fetch (:chr vc) (:start vc) (:end vc))
                      vc-type (if-not (empty? other-vcs) :total :unique)]
                  (assoc coll vc-type (inc (get coll vc-type)))))
              {:total 0 :unique 0}
              (parse-vcf file1)))))

(defn top-level-metrics
  "Provide one-line summary of similarity metrics for a VCF comparison."
  [x]
  (letfn [(passes-filter? [vc]
            (= (count (:filters vc)) 0))
          (nonref-passes-filter? [vc]
            (and (passes-filter? vc)
                 (every? #(contains? #{"HET" "HOM_VAR"} (:type %)) (:genotypes vc))))
          (vrn-type-passes-filter [vrn-type]
            (fn [vc]
              (and (passes-filter? vc)
                   (contains? vrn-type (:type vc)))))
          (count-variants [f check?]
            (count (filter check? (parse-vcf f))))
          (all-vrn-counts [fname]
            {:total (count-variants fname passes-filter?)
             :snp (count-variants fname (vrn-type-passes-filter #{"SNP"}))
             :indel (count-variants fname (vrn-type-passes-filter #{"INDEL"}))})]
    (ordered-map
     :sample (:sample x)
     :call1 (-> x :c1 :name)
     :call2 (-> x :c2 :name)
     :genotype_concordance (-> x :metrics :percent_overall_genotype_concordance)
     :nonref_discrepency (-> x :metrics :percent_non_reference_discrepancy_rate)
     :nonref_sensitivity (-> x :metrics :percent_non_reference_sensitivity)
     :concordant (all-vrn-counts (first (:c-files x)))
     :nonref_concordant (count-variants (first (:c-files x)) nonref-passes-filter?)
     :discordant1 (all-vrn-counts (second (:c-files x)))
     :discordant2 (all-vrn-counts (nth (:c-files x) 2))
     :discordant_both (apply discordance-metrics (rest (:c-files x))))))

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
                              [:discordant1 :snp] (str "Discordant SNPs: " (:call1 metrics))
                              [:discordant1 :indel] (str "Discordant indels: " (:call1 metrics))
                              [:discordant2 :snp] (str "Discordant SNPs: " (:call2 metrics))
                              [:discordant2 :indel] (str "Discordant indels: " (:call2 metrics)))]
    (letfn [(get-value-and-score [[k metric]]
              (let [val (if (coll? k) (get-in metrics k) (get metrics k))]
                {:metric metric :value val :score (calc-score k val)}))]
      (.write wrtr (str (doric/table [:metric :value]
                                     (map get-value-and-score to-write))
                        "\n")))))
