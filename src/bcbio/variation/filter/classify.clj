(ns bcbio.variation.filter.classify
  "Provide classification based filtering for variants."
  (:import [org.broadinstitute.sting.utils.variantcontext VariantContextBuilder]
           [org.broadinstitute.sting.utils.codecs.vcf VCFHeader VCFInfoHeaderLine
            VCFHeaderLineType VCFFilterHeaderLine])
  (:use [ordered.set :only [ordered-set]]
        [clj-ml.utils :only [serialize-to-file deserialize-from-file]]
        [clj-ml.data :only [make-dataset dataset-set-class make-instance]]
        [clj-ml.classifiers :only [make-classifier classifier-train
                                   classifier-evaluate classifier-classify]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template
                                               get-vcf-source]])
  (:require [clojure.string :as string]
            [incanter.stats :as stats]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]))

;; ## Normalized attribute access

(defmulti get-vc-attr
  "Generalized retrieval of attributes from variant with a single genotype."
  (fn [vc attr] attr))

(defmethod get-vc-attr "AD"
  [vc attr]
  "AD: Allelic depth for ref and alt alleles. Converted to percent
   deviation from expected for haploid/diploid calls.
   Also calculates allele depth from AO and DP used by FreeBayes.
   AO is the count of the alternative allele."
  (letfn [(calc-expected [g ref-count allele-count]
            {:pre [(not (neg? ref-count))]}
            (when-let [e-pct (get {"HOM_VAR" 1.0 "HET" 0.5 "HOM_REF" 0.0} (:type g))]
              (Math/abs (- e-pct (/ allele-count (+ allele-count ref-count))))))
          (from-ad [g]
            (let [ads (map #(Float/parseFloat %) (string/split (get-in g [:attributes attr]) #","))
                  alleles (cons (:ref-allele vc) (:alt-alleles vc))
                  ref-count (first ads)
                  allele-count (apply + (map #(nth ads (.indexOf alleles %)) (set (:alleles g))))]
              (calc-expected g ref-count allele-count)))
          (from-ao [g]
            (let [alt-count (Float/parseFloat (get-in g [:attributes "AO"]))
                  total-count (Float/parseFloat (get-in g [:attributes "DP"]))]
              (calc-expected g (- total-count alt-count) alt-count)))]
    (let [g (-> vc :genotypes first)]
      (cond
       (get-in g [:attributes attr]) (from-ad g)
       (get-in g [:attributes "AO"]) (from-ao g)
       :else (Exception. (str "AD not found in attributes" (:attributes g)))))))

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
  {:pre [(= 1 (:num-samples vc))
         (contains? #{1 2} (-> vc :genotypes first :alleles count))]}
  (zipmap attrs (map (partial get-vc-attr vc) attrs)))

(defn get-vc-attr-ranges
  "Retrieve quantile ranges of attributes for min/max normalization."
  [attrs in-vcf ref]
  (letfn [(get-quartiles [[k v]]
            [k (stats/quantile v :probs [0.05 0.95])])]
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
            (let [safe-maxv (if (= minv maxv) (inc maxv) maxv)
                  trunc-score-max (if (< x safe-maxv) x safe-maxv)
                  trunc-score (if (> trunc-score-max minv) trunc-score-max minv)]
              (/ (- trunc-score minv) (- safe-maxv minv))))
          (min-max-norm-ranges [mm-ranges [k v]]
            [k (min-max-norm v (get mm-ranges k))])]
    (let [mm-ranges (get-vc-attr-ranges attrs in-vcf ref)]
      (fn [vc]
        (->> (get-vc-attrs vc attrs)
             (map (partial min-max-norm-ranges mm-ranges))
             (into {}))))))

;; ## Linear classifier

(defn- get-vc-inputs
  [attrs normalizer group vc]
  (let [n-vals (normalizer vc)]
    (conj (vec (map #(get n-vals %) attrs)) group)))

(defn- get-train-inputs
  "Retrieve normalized training inputs from VCF file."
  [group in-vcf attrs normalizer ref]
  (with-open [vcf-s (get-vcf-source in-vcf ref)]
    (doall (map (partial get-vc-inputs attrs normalizer group)
                (parse-vcf vcf-s)))))

(defn- train-vcf-classifier
  "Do the work of training a variant classifier."
  [attrs base-vcf true-vcf false-vcf ref]
  (let [normalizer (get-vc-attrs-normalized attrs base-vcf ref)]
    (-> (make-classifier :regression :linear)
        (classifier-train (make-dataset "ds" (conj attrs :c)
                                        (concat (get-train-inputs 1 true-vcf attrs normalizer ref)
                                                (get-train-inputs 0 false-vcf attrs normalizer ref))
                                        {:class :c})))))

(defn build-vcf-classifier
  "Provide a variant classifier based on provided attributes and true/false examples."
  [attrs base-vcf true-vcf false-vcf ref]
  (let [out-file (str (itx/file-root base-vcf) "-classifier.bin")]
    (if-not (itx/needs-run? out-file)
      (deserialize-from-file out-file)
      (let [classifier (train-vcf-classifier attrs base-vcf true-vcf false-vcf ref)]
        (serialize-to-file classifier out-file)
        classifier))))

(defn- add-cfilter-header
  "Add details on the filtering to the VCF file header."
  [attrs]
  (fn [header]
    (let [desc (str "Classification score based on true/false positives for: "
                    (string/join ", " attrs))
          new #{(VCFInfoHeaderLine. "CSCORE" 1 VCFHeaderLineType/Float desc)
                (VCFFilterHeaderLine. "CScoreFilter" "Based on classifcation CSCORE")}]
      (VCFHeader. (apply ordered-set (concat (.getMetaData header) new))
                  (.getGenotypeSamples header)))))

(defn- filter-vc
  "Update a variant context with filter information from classifier."
  [classifier normalizer config vc]
  (let [attrs (vec (:classifiers config))
        score (-> (make-dataset "ds" (conj attrs :c)
                                 [(get-vc-inputs attrs normalizer -1 vc)]
                                 {:class :c})
                  (make-instance (assoc (normalizer vc) :c -1))
                  (#(classifier-classify classifier %)))]
    (-> (VariantContextBuilder. (:vc vc))
        (.attributes (assoc (:attributes vc) "CSCORE" score))
        (.filters (when (< score (:min-cscore config))
                    #{"CScoreFilter"}))
        .make)))

(defn filter-vcf-w-classifier
  "Filter an input VCF file using a trained classifer on true/false variants."
  [base-vcf true-vcf false-vcf ref config]
  (let [out-file (itx/add-file-part base-vcf "cfilter")
        c (build-vcf-classifier (:classifiers config) base-vcf
                                true-vcf false-vcf ref)
        normalizer (get-vc-attrs-normalized (:classifiers config) base-vcf ref)]
    (when (itx/needs-run? out-file)
      (with-open [vcf-s (get-vcf-source base-vcf ref)]
        (write-vcf-w-template base-vcf {:out out-file}
                              (map (partial filter-vc c normalizer config)
                                   (parse-vcf vcf-s))
                              ref :header-update-fn (add-cfilter-header (:classifiers config)))))
    out-file))

(defn pipeline-classify-filter
  "Fit VCF classification-based filtering into analysis pipeline."
  [in-vcf train-info ref config]
  (letfn [(get-train-vcf [type]
            (-> (filter #(= type (:name %)) train-info)
                       first
                       :file))]
    (filter-vcf-w-classifier in-vcf (get-train-vcf "concordant")
                             (get-train-vcf "discordant") ref config)))
