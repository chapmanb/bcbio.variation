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
        [bcbio.variation.haploid :only [get-likelihoods]]
        [bcbio.variation.metrics :only [to-float]]
        [bcbio.variation.filter.intervals :only [pipeline-combine-intervals]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template
                                               get-vcf-iterator has-variants?
                                               get-vcf-retriever variants-in-region]])
  (:require [clojure.string :as string]
            [incanter.stats :as stats]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.variation.index.gemini :as gemini]))

;; ## Normalized attribute access

(defmulti get-vc-attr
  "Generalized retrieval of attributes from variant with a single genotype."
  (fn [vc attr retrievers]
    (if (contains? gemini/gemini-metrics attr)
      :gemini
      attr)))

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
            (let [alt-count (Float/parseFloat (get-in g [:attributes "AO"]))
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
     (zero? (count pls)) 0.0
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
    (let [x (getter vc attr)]
      (if-not (nil? x) x
              (cond
               (.startsWith attr "gms") 100.0
               :else nil)))))

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
    (with-open [vcf-iter (get-vcf-iterator in-vcf ref)]
      (->> (reduce (fn [coll vc]
                    (reduce (fn [icoll [k v]]
                              (assoc icoll k (cons v (get icoll k))))
                            coll (get-vc-attrs vc attrs retrievers)))
                  (zipmap attrs (repeat [])) (parse-vcf vcf-iter))
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

;; ## Linear classifier

(defn- get-vc-inputs
  [attrs normalizer group vc]
  (let [n-vals (normalizer vc)]
    (conj (vec (map #(get n-vals %) attrs)) group)))

(def ^{:private true
       :doc "Available specialized classifier groups."}
  classifier-types [:snp :complex])

(defn- get-classifier-type
  "Map variant types to specialized classifiers."
  [vc]
  (case (:type vc)
    "SNP" :snp
    :complex))

(defn- get-train-inputs
  "Retrieve normalized training inputs from VCF file."
  [group in-vcf ctype attrs normalizer ref]
  (with-open [vcf-iter (get-vcf-iterator in-vcf ref)]
    (->> (parse-vcf vcf-iter)
         (filter #(= ctype (get-classifier-type %)))
         (map (partial get-vc-inputs attrs normalizer group))
         doall)))

(defn- train-vcf-classifier
  "Do the work of training a variant classifier."
  [ctype attrs base-vcf true-vcf false-vcf ref config]
  (let [normalizer (partial get-vc-attrs-normalized attrs base-vcf ref config)
        ds (make-dataset "ds" (conj (vec attrs) {:c [:a :b]})
                         (concat (get-train-inputs :a true-vcf ctype attrs
                                                   (normalizer true-vcf)
                                                   ref)
                                 (get-train-inputs :b false-vcf ctype attrs
                                                   (normalizer false-vcf)
                                                   ref))
                      {:class :c})
        c (classifier-train (make-classifier :support-vector-machine :smo) ds)]
    ;(println "Evaluate" (classifier-evaluate c :dataset ds ds))
    c))

(defn- build-vcf-classifiers
  "Provide a variant classifier based on provided attributes and true/false examples."
  [attrs base-vcf true-vcf false-vcf ref config]
  (letfn [(build-vcf-classifier [ctype]
            (let [out-file (format "%s-%s-classifier.bin" (itx/file-root base-vcf) (name ctype))]
              (if-not (itx/needs-run? out-file)
                (deserialize-from-file out-file)
                (let [classifier (train-vcf-classifier ctype attrs base-vcf true-vcf false-vcf
                                                       ref config)]
                  (serialize-to-file classifier out-file)
                  classifier))))]
    (zipmap classifier-types (map build-vcf-classifier classifier-types))))

(defn- add-cfilter-header
  "Add details on the filtering to the VCF file header."
  [attrs]
  (fn [_ header]
    (let [desc (str "Classification score based on true/false positives for: "
                    (string/join ", " attrs))
          new #{(VCFInfoHeaderLine. "CSCORE" 1 VCFHeaderLineType/Float desc)
                (VCFFilterHeaderLine. "CScoreFilter" "Based on classifcation CSCORE")}]
      (VCFHeader. (apply ordered-set (concat (.getMetaDataInInputOrder header) new))
                  (.getGenotypeSamples header)))))

(defn- filter-vc
  "Update a variant context with filter information from classifier."
  [classifiers normalizer trusted-get config vc]
  (let [attrs (vec (:classifiers config))
        score (-> (make-dataset "ds" (conj attrs :c)
                                 [(get-vc-inputs attrs normalizer -1 vc)]
                                 {:class :c})
                  (make-instance (assoc (normalizer vc) :c -1))
                  (#(classifier-classify (get classifiers (get-classifier-type vc)) %)))]
    (-> (VariantContextBuilder. (:vc vc))
        (.attributes (assoc (:attributes vc) "CSCORE" score))
        (.filters (when (and (not (has-variants? trusted-get
                                                 (:chr vc) (:start vc) (:end vc)
                                                 (:ref-allele vc) (:alt-alleles vc)))
                             (< score (get config :min-cscore 0.5)))
                    #{"CScoreFilter"}))
        .make)))

(defn- get-target-variants
  "Retrieve variants from original file based on variants in target file."
  [orig-file target-file ref-file ext]
  (letfn [(get-orig-variants [retriever vc]
            (->> (variants-in-region retriever (:chr vc) (:start vc) (:end vc))
                 (filter #(= (:start %) (:start vc)))
                 (map :vc)))]
    (let [out-file (itx/add-file-part orig-file ext)]
      (when (itx/needs-run? out-file)
        (with-open [vcf-iter (get-vcf-iterator target-file ref-file)
                    retriever (get-vcf-retriever ref-file orig-file)]
          (write-vcf-w-template orig-file {:out out-file}
                                (mapcat (partial get-orig-variants retriever)
                                        (parse-vcf vcf-iter))
                                ref-file)))
      out-file)))

(defn filter-vcf-w-classifier
  "Filter an input VCF file using a trained classifier on true/false variants."
  [base-vcf orig-true-vcf orig-false-vcf trusted-vcf ref config & {:keys [train-w-base?]}]
  (let [out-file (itx/add-file-part base-vcf "cfilter")
        true-vcf (if train-w-base?
                   (get-target-variants base-vcf orig-true-vcf ref "tps")
                   orig-true-vcf)
        false-vcf (if train-w-base?
                    (get-target-variants base-vcf orig-false-vcf ref "fps")
                    orig-false-vcf)
        cs (build-vcf-classifiers (:classifiers config) base-vcf
                                  true-vcf false-vcf ref config)
        normalizer (get-vc-attrs-normalized (:classifiers config) base-vcf ref config base-vcf)]
    (when (itx/needs-run? out-file)
      (println "Filter VCF with" (str cs))
      (with-open [vcf-iter (get-vcf-iterator base-vcf ref)
                  trusted-get (get-vcf-retriever ref trusted-vcf)]
        (write-vcf-w-template base-vcf {:out out-file}
                              (map (partial filter-vc cs normalizer trusted-get config)
                                   (parse-vcf vcf-iter))
                              ref :header-update-fn (add-cfilter-header (:classifiers config)))))
    out-file))

(defn pipeline-classify-filter
  "Fit VCF classification-based filtering into analysis pipeline."
  [in-vcf train-info call exp params config]
  (letfn [(get-train-vcf [type]
            (-> (filter #(= type (:name %)) train-info)
                       first
                       :file))]
    (pipeline-combine-intervals exp config)
    (filter-vcf-w-classifier in-vcf (get-train-vcf "concordant")
                             (get-train-vcf "discordant")
                             (get-train-vcf "trusted")
                             (:ref exp) params
                             :train-w-base? (get call :recall false))))
