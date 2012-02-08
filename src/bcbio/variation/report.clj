;; Parse and provide detailed information from GATKReport output

(ns bcbio.variation.report
  (:import [org.broadinstitute.sting.gatk.report GATKReport])
  (:use [ordered.map :only [ordered-map]]
        [clojure.math.combinatorics :only [cartesian-product]]
        [bcbio.variation.variantcontext :only [parse-vcf]])
  (:require [doric.core :as doric]))

(defn concordance-report-metrics [sample in-file]
  "Retrieve high level concordance metrics from GATK VariantEval report."
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

(defn top-level-metrics [x]
  "Provide one-line summary of similarity metrics for a VCF comparison."
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
     :discordant2 (all-vrn-counts (nth (:c-files x) 2)))))

(defn calc-accuracy [metrics]
  "Calculate an overall accuracy score from input metrics.
  The accuracy logic is:
  (#correctly aligned bases / (#correctly aligned bases +
                               1*(simple substitutions and indels) +
                               2*(larger errors)))."
  (letfn [(get-penalty [[error-type call-type]]
            (case call-type
              :snp 1
              :indel 2))]
    (let [error-items (cartesian-product [:discordant :phasing-error] [:snp :indel])
          error-score (apply + (map #(* (get-in metrics %) (get-penalty %)) error-items))]
      (float
       (/ (:total-bases metrics)
          (+ (:total-bases metrics) error-score))))))

(defn write-scoring-table [metrics wrtr]
  "Summary table of high level variables and scoring metrics for comparison."
  (let [to-write (ordered-map :accuracy "Overall accuracy score"
                              :total-bases "Total bases compared"
                              [:discordant :snp] "Discordant SNPs"
                              [:discordant :indel] "Discordant indels"
                              [:phasing-error :snp] "Phasing Error SNPs"
                              [:phasing-error :indel] "Phasing Error indels"
                              :haplotype-blocks "Phased haplotype blocks"
                              :nonmatch-het-alt "Non-matching heterozygous alternative alleles")
        s-metrics (assoc metrics :accuracy (calc-accuracy metrics))]
    (letfn [(prep-row [[k x]]
              {:metric x
               :value (if (coll? k) (get-in s-metrics k) (get s-metrics k))})]
      (.write wrtr (str (doric/table [:metric :value] (map prep-row to-write))
                        "\n")))))

(defn calc-score [type val]
  "Calculate scoring for input metric types"
  (if (keyword? type) ""
      (let [sign (if (= :concordant (first type)) 1 -1)]
        (case (second type)
              :snp (-> 1 (* sign) (* val))
              :indel (-> 2 (* sign) (* val))
              (throw (Exception. (str "Unexpected variant type" type)))))))

(defn write-concordance-metrics [metrics wrtr]
  "Summary table of metrics for assessing the score of a variant comparison."
  (let [to-write (ordered-map :genotype_concordance "Overall genotype concordance"
                              :nonref_discrepency "Non-reference discrepancy rate"
                              :nonref_sensitivity "Non-reference sensitivity"
                              :nonref_concordant "Non-reference concordant count"
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
