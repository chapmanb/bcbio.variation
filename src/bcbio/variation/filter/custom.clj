(ns bcbio.variation.filter.custom
  "Provide custom filtering tied to specific variant callers."
  (:require [clojure.walk :refer [keywordize-keys]]
            [bcbio.variation.filter.attr :as attr]
            [bcbio.variation.variantcontext :as gvc]))

(defn- freebayes-filter?
  "Define filtered FreeBayes calls based on call type, depth and two measures of strand bias.
   Uses approach determined from manual inspection of true/false calls on
   NA12878 compared to Genome in a Bottle Reference Genomes.
   Filters differently for homozygote and heterozygote variants based on
   depth, quality and two strand bias metrics:
   - AD: deviation from expected bias.
   - QR_QA: Sequence quality relationships between reference and alternative alleles."
  [vc]
  {:pre [(= 1 (count (:genotypes vc)))]}
  (let [attrs (keywordize-keys (attr/get-vc-attrs vc ["DP" "QUAL" "AD"] {}))]
    (case (-> vc :genotypes first :type)
      "HET" (or (< (:DP attrs) 4)
                (and (< (:DP attrs) 13)
                     (< (:QUAL attrs) 20)
                     (> (:AD attrs) 0.1)))
      "HOM_VAR" (or (and (< (:DP attrs) 4) (< (:QUAL attrs) 50))
                    (and (< (:DP attrs) 13) (> (:AD attrs) 0.1)))
      false)))

(defn freebayes-filter
  "Custom filtering of FreeBayes handling depth, strand bias and quality."
  [vcf-file ref-file]
  (gvc/write-vcf-from-filter vcf-file ref-file "filter" "FreeBayesCustom"
                             "bcbio.variation FreeBayes filter: DP, QUAL and strand bias"
                             (complement freebayes-filter?)))

(def ^{:doc "Available custom filters" :private true} cfilters
  {:freebayes #'freebayes-filter})

(defn -main [& args]
  (if-let [filter-fn (get cfilters (keyword (first args)))]
    (apply filter-fn (rest args))
    (do
      (println "Incorrect filter name. Available filters:")
      (doseq [k (sort (keys cfilters))]
        (println (format " %s  %s" (name k) (:doc (meta (get cfilters k)))))))))
