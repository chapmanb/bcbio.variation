(ns bcbio.variation.filter.attr
  "Provide generalized access to variant attributes, handling retrieval
   from multiple sources (VCF INFO file, VCF FORMAT field, Gemini)."
  (:use [bcbio.variation.haploid :only [get-likelihoods]]
        [bcbio.variation.metrics :only [to-float]]
        )
  (:require [clojure.string :as string]
            [incanter.stats :as stats]
            [bcbio.variation.variantcontext :as gvc]
            [bcbio.variation.index.gemini :as gemini]))

(defmulti get-vc-attr
  "Generalized retrieval of attributes from variant with a single genotype."
  (let [gemini-ids (set (map :id (gemini/available-metrics nil :include-noviz? true)))]
    (fn [vc attr retrievers]
      (if (contains? gemini-ids attr)
        :gemini
        attr))))

(defmethod get-vc-attr "AD"
  ^{:doc "AD: Allelic depth for ref and alt alleles. Converted to percent
          deviation from expected for haploid/diploid calls.
          Also calculates allele depth from AO and DP used by FreeBayes.
          AO is the count of the alternative allele."}
  [vc attr _]
  {:pre [(= 1 (:num-samples vc))
         (contains? #{1 2} (-> vc :genotypes first :alleles count))]}
  (letfn [(calc-expected [g ref-count allele-count]
            {:pre [(not (neg? ref-count))]}
            (when (or (pos? ref-count) (pos? allele-count))
              (when-let [e-pct (get {"HOM_VAR" 1.0 "HET" 0.5 "HOM_REF" 0.0} (:type g))]
                (Math/abs (- e-pct (/ allele-count (+ allele-count ref-count)))))))
          (from-ad [g]
            (let [ads (map float (get-in g [:attributes attr]))
                  alleles (cons (:ref-allele vc) (:alt-alleles vc))
                  ref-count (first ads)
                  allele-count (apply + (map #(nth ads (.indexOf alleles %)) (set (:alleles g))))]
              (calc-expected g ref-count allele-count)))
          (from-ao [g]
            (let [alt-count (apply max (map #(Float/parseFloat %) (string/split (get-in g [:attributes "AO"]) #",")))
                  total-count (float (get-in g [:attributes "DP"]))]
              (calc-expected g (- total-count alt-count) alt-count)))]
    (let [g (-> vc :genotypes first)]
      (cond
       (get-in g [:attributes "AO"]) (from-ao g)
       (seq (get-in g [:attributes attr])) (from-ad g)
       :else nil
       ;; (println (format "AD not found in attributes %s %s %s"
       ;;                  (:attributes g) (:chr vc) (:start vc)))
       ))))

(defmethod get-vc-attr [:format "AD"]
  [vc attr retrievers]
  (get-vc-attr vc "AD" retrievers))

(defmethod get-vc-attr "PL"
  ^{:doc "Provide likelihood ratios for genotype compared to next most likely call.
          For haploid calls, get the homozygous reference or homozygous variant
          likelihood. For diploid calls, get the largest alternative value."}
  [vc attr _]
  {:pre [(= 1 (:num-samples vc))
         (contains? #{1 2} (-> vc :genotypes first :alleles count))]}
  (let [g (-> vc :genotypes first)
        pls (dissoc (get-likelihoods (:genotype g) :no-convert true)
                    (:type g))]
    (cond
     (zero? (count pls)) nil
     (> (count (:alleles g)) 1) (apply max (vals pls))
     (= (:type g) "HOM_VAR") (get pls "HOM_REF")
     (= (:type g) "HOM_REF") (get pls "HOM_VAR"))))

(defmethod get-vc-attr "QUAL"
  [vc attr _]
  (:qual vc))

(defmethod get-vc-attr [:format "DP"]
  ^{:doc "Retrieve depth from Genotype FORMAT metrics.
          Handles custom cases like cortex_var with alternative
          depth attributes."}
  [vc attr _]
  {:pre [(= 1 (:num-samples vc))]}
  (let [g-attrs (-> vc :genotypes first :attributes)]
    (cond
     (contains? g-attrs "COV") (int (apply + (map to-float (string/split (get g-attrs "COV") #","))))
     (contains? g-attrs "DP") (get g-attrs "DP")
     (contains? g-attrs "AD") (apply + (get g-attrs "AD"))
     :else nil)))

(defmethod get-vc-attr :gemini
  ^{:doc "Retrieve attribute information from associated Gemini index."}
  [vc attr retrievers]
  (when-let [getter (:gemini retrievers)]
    (getter vc attr)))

(defmethod get-vc-attr :default
  [vc attr _]
  (let [x (get-in vc [:attributes attr])]
    (when-not (nil? x)
      (try (Float/parseFloat x)
           (catch java.lang.NumberFormatException _ x)))))

(defn get-vc-attrs
  "Retrieve attributes from variants independent of location."
  [vc attrs retrievers]
  (zipmap attrs (map #(get-vc-attr vc % retrievers) attrs)))

(defn get-vc-attr-ranges
  "Retrieve quantile ranges of attributes for min/max normalization."
  [attrs in-vcf ref retrievers]
  (letfn [(get-quartiles [[k v]]
            [k (stats/quantile v :probs [0.05 0.95])])]
    (with-open [vcf-iter (gvc/get-vcf-iterator in-vcf ref)]
      (->> (reduce (fn [coll vc]
                    (reduce (fn [icoll [k v]]
                              (assoc icoll k (cons v (get icoll k))))
                            coll (get-vc-attrs vc attrs retrievers)))
                  (zipmap attrs (repeat [])) (gvc/parse-vcf vcf-iter))
           (map get-quartiles)
           (into {})))))

(defn- get-external-retrievers
  [in-file ref-file]
  {:gemini (gemini/vc-attr-retriever in-file ref-file)})

(defmulti get-vc-attrs-normalized
  "Normalized attributes for each variant context in an input file.
   Passed two input VCFs:
    - in-vcf -- provides full range of inputs for classification and
      used for building normalization ranges.
    - work-vcf -- file for attribute retrieval, used to setup variable
      retrieval from external sources like Gemini"
  (fn [_ _ _ config _] (keyword (get config :normalize "default"))))

;; Min-max normalization
(defmethod get-vc-attrs-normalized :minmax
  [attrs in-vcf ref config work-vcf]
  (letfn [(min-max-norm [x [minv maxv]]
            (let [safe-maxv (if (= minv maxv) (inc maxv) maxv)
                  trunc-score-max (if (< x safe-maxv) x safe-maxv)
                  trunc-score (if (> trunc-score-max minv) trunc-score-max minv)]
              (/ (- trunc-score minv) (- safe-maxv minv))))
          (min-max-norm-ranges [mm-ranges [k v]]
            [k (min-max-norm v (get mm-ranges k))])]
    (let [retrievers (get-external-retrievers in-vcf ref)
          mm-ranges (get-vc-attr-ranges attrs in-vcf ref retrievers)
          work-retrievers (get-external-retrievers work-vcf ref)]
      (fn [vc]
        (->> (get-vc-attrs vc attrs work-retrievers)
             (map (partial min-max-norm-ranges mm-ranges))
             (into {}))))))

;; No normalization
(defmethod get-vc-attrs-normalized :default
  [attrs _ ref config work-vcf]
  (let [retrievers (get-external-retrievers work-vcf ref)]
    (fn [vc]
      (into {} (get-vc-attrs vc attrs retrievers)))))

(defn prep-vc-attr-retriever
  "Provide easy lookup of attributes from multiple input sources"
  [in-file ref-file]
  (let [retrievers (get-external-retrievers in-file ref-file)]
    (fn [attrs vc]
      (into {} (get-vc-attrs vc attrs retrievers)))))
