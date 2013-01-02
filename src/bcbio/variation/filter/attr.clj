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
                  ref-count (first ads)
                  allele-count (apply + (rest ads))]
              (calc-expected g ref-count allele-count)))
          (from-ao [g]
            (let [alt-count (apply + (map #(Float/parseFloat %)
                                          (string/split (get-in g [:attributes "AO"]) #",")))
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
  ^{:doc "Provide likelihood confidence for the called genotype.
          For reference calls, retrieve the likelihood of the most likely
          variant (least negative). For variant calls, retrieve
          the reference likelihood."}
  [vc attr _]
  {:pre [(= 1 (:num-samples vc))
         (contains? #{1 2} (-> vc :genotypes first :alleles count))]}
  (let [g (-> vc :genotypes first)
        pls (dissoc (get-likelihoods (:genotype g) :no-convert true)
                    (:type g))]
    (cond
     (zero? (count pls)) nil
     (= (:type g) "HOM_REF") (apply max (vals pls))
     :else (get pls "HOM_REF"))))

(defmethod get-vc-attr "PL-ratio"
  ^{:doc "Calculate ratio of reference likelihood call to alternative variant calls.
          This helps measure whether a call is increasingly likely to be reference
          compared with variant choices."}
  [vc attr _]
  {:pre [(= 1 (:num-samples vc))]}
  (let [g (-> vc :genotypes first)
        pls (dissoc (get-likelihoods (:genotype g) :no-convert true)
                    (:type g))]
    (when-not (zero? (count pls))
      (/ (get pls "HOM_REF")
         (apply min (cons -1.0 (-> pls (dissoc "HOM_REF") vals)))))))

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
     (contains? g-attrs "DP") (to-float (get g-attrs "DP"))
     (contains? g-attrs "AD") (to-float (apply + (get g-attrs "AD")))
     (contains? g-attrs "COV") (int (apply + (map to-float (string/split (get g-attrs "COV") #","))))
     :else nil)))

(defmethod get-vc-attr "DP"
  ^{:doc "Retrieve depth for an allele, first trying genotype information
          then falling back on information in INFO column."}
  [vc attr rets]
  (let [gt-dp (when (= 1 (:num-samples vc))
                (get-vc-attr vc [:format attr] rets))
        info-dp (let [x (get-in vc [:attributes attr])]
                  (when-not (nil? x)
                    (try (Float/parseFloat x)
                         (catch java.lang.NumberFormatException _ x))))
        ready-dp (if (nil? gt-dp) info-dp gt-dp)]
    (if (nil? gt-dp) info-dp gt-dp)))

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
            [k (stats/quantile (remove nil? v) :probs [0.05 0.95])])]
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
  (fn [_ _ _ config] (keyword (get config :normalize "default"))))

(defmethod get-vc-attrs-normalized :minmax
  ^{:doc "Minimum/maximum normalization to a 0-1 scale using quartiles."}
  [attrs in-vcf ref config]
  (letfn [(min-max-norm [x [minv maxv]]
            (let [safe-maxv (if (= minv maxv) (inc maxv) maxv)
                  trunc-score-max (if (< x safe-maxv) x safe-maxv)
                  trunc-score (if (> trunc-score-max minv) trunc-score-max minv)]
              (/ (- trunc-score minv) (- safe-maxv minv))))
          (min-max-norm-ranges [mm-ranges [k v]]
            [k (when-not (nil? v)
                 (min-max-norm v (get mm-ranges k)))])]
    (let [retrievers (get-external-retrievers in-vcf ref)
          mm-ranges (get-vc-attr-ranges attrs in-vcf ref retrievers)]
      (fn [work-vcf]
        (let [work-retrievers (get-external-retrievers work-vcf ref)]
          (fn [vc]
            (->> (get-vc-attrs vc attrs work-retrievers)
                 (map (partial min-max-norm-ranges mm-ranges))
                 (into {}))))))))

(defmethod get-vc-attrs-normalized :log
  ^{:doc "Log normalization of specified input variables."}
  [attrs in-vcf ref config]
  (let [base-fn (get-vc-attrs-normalized attrs in-vcf ref (assoc config :normalize :default))
        need-log-vars (set (:log-attrs config))]
    (fn [work-vcf]
      (let [inner-fn (base-fn work-vcf)]
        (fn [vc]
          (reduce (fn [coll [k v]]
                    (assoc coll k
                           (if (contains? need-log-vars k) (Math/log v) v)))
                  {} (inner-fn vc)))))))

(defmethod get-vc-attrs-normalized :default
  ^{:doc "Attribute access without normalization."}
  [attrs _ ref config]
  (fn [work-vcf]
    (let [retrievers (get-external-retrievers work-vcf ref)]
      (fn [vc]
        (into {} (get-vc-attrs vc attrs retrievers))))))

(defn prep-vc-attr-retriever
  "Provide easy lookup of attributes from multiple input sources"
  [in-file ref-file]
  (let [retrievers (get-external-retrievers in-file ref-file)]
    (fn [attrs vc]
      (into {} (get-vc-attrs vc attrs retrievers)))))
