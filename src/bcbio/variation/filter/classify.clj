(ns bcbio.variation.filter.classify
  "Provide classification based filtering for variants."
  (:use [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template
                                               get-vcf-source]])
  (:require [clojure.string :as string]
            [incanter.stats :as stats]
            [bcbio.run.itx :as itx]))

;; ## Normalized attribute access

(defmulti get-vc-attr
  "Generalized retrieval of attributes from variant with a single genotype."
  (fn [vc attr] attr))

(defmethod get-vc-attr "AD"
  [vc attr]
  "AD: Allelic depth for ref and alt alleles. Converted to percent
   deviation from expected for haploid/diploid calls."
  (let [g (-> vc :genotypes first)
        ads (map #(Integer/parseInt %) (string/split (get-in g [:attributes attr]) #","))
        alleles (cons (:ref-allele vc) (:alt-alleles vc))
        ref-count (first ads)
        allele-count (apply + (map #(nth ads (.indexOf alleles %)) (set (:alleles g))))]
    (when-let [e-pct (get {"HOM_VAR" 1.0 "HET" 0.5 "HOM_REF" 0.0} (:type g))]
      (Math/abs (- e-pct (/ allele-count (+ allele-count ref-count)))))))

(defmethod get-vc-attr "QUAL"
  [vc attr]
  (:qual vc))

(defmethod get-vc-attr :default
  [vc attr]
  (let [x (get-in vc [:attributes attr])]
    (try (Float/parseFloat x)
         (catch java.lang.NumberFormatException _ x))))

(defn get-vc-attrs
  "Retrieve attributes from variants independent of location."
  [vc attrs]
  {:pre [(= 1 (count (:genotypes vc)))
         (contains? #{1 2} (-> vc :genotypes first :alleles count))]}
  (zipmap attrs (map (partial get-vc-attr vc) attrs)))

(defn get-vc-attr-ranges
  "Retrieve first/third quantile ranges of attributes for min/max normalization."
  [attrs in-vcf ref]
  (letfn [(get-quartiles [[k v]]
            [k (stats/quantile v :probs [0.25 0.75])])]
    (with-open [vcf-s (get-vcf-source in-vcf ref)]
      (->> (reduce (fn [coll vc]
                    (reduce (fn [icoll [k v]]
                              (assoc icoll k (cons v (get icoll k))))
                            coll (get-vc-attrs vc attrs)))
                  (zipmap attrs (repeat [])) (parse-vcf vcf-s))
           (map get-quartiles)
           (into {})))))

(defn get-vc-attrs-normalized
  "Min-max Normalized attributes for each variant context in an input file."
  [attrs in-vcf ref]
  (letfn [(min-max-norm [x [minv maxv]]
            (let [trunc-score-max (if (< x maxv) x maxv)
                  trunc-score (if (> trunc-score-max minv) trunc-score-max minv)]
              (/ (- trunc-score minv) (- maxv minv))))
          (min-max-norm-ranges [mm-ranges [k v]]
            [k (min-max-norm v (get mm-ranges k))])]
    (let [mm-ranges (get-vc-attr-ranges attrs in-vcf ref)]
      (fn [vc]
        (->> (get-vc-attrs vc attrs)
             (map (partial min-max-norm-ranges mm-ranges))
             (into {}))))))
