(ns bcbio.variation.filter.classify
  "Provide classification based filtering for variants."
  (:import [org.broadinstitute.sting.utils.variantcontext VariantContextBuilder]
           [org.broadinstitute.sting.utils.codecs.vcf VCFHeader VCFInfoHeaderLine
            VCFHeaderLineType VCFFilterHeaderLine])
  (:use [ordered.set :only [ordered-set]]
        [clojure.math.combinatorics :only [cartesian-product]]
        [clj-ml.utils :only [serialize-to-file deserialize-from-file]]
        [clj-ml.data :only [make-dataset dataset-set-class make-instance]]
        [clj-ml.classifiers :only [make-classifier classifier-train
                                   classifier-evaluate classifier-classify]]
        [bcbio.variation.filter.attr :only [get-vc-attrs-normalized prep-vc-attr-retriever]]
        [bcbio.variation.filter.intervals :only [pipeline-combine-intervals]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template
                                               get-vcf-iterator has-variants?
                                               get-vcf-retriever variants-in-region]])
  (:require [clojure.string :as string]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]))

;; ## Split variants for classification

(defn- classifier-types
  "Define splitting of classifiers based on variant characteristics."
  []
  (let [variant-types [:snp :complex]
        repeats [true false]]
    (map (fn [[vtype rpt]] {:variant-type vtype
                            :repetitive rpt})
         (cartesian-product variant-types repeats))))

(defn- ctype-to-str
  "Convert a classifier types into a string name for output files."
  [x]
  (str (name (:variant-type x))
       (if (:repetitive x) "rpt" "std")))

(defn- get-classifier-type
  "Map variant types to specialized classifiers."
  [vc attr-get]
  (let [attrs (attr-get ["rmsk"] vc)]
    {:variant-type (case (:type vc)
                     "SNP" :snp
                     :complex)
     :repetitive (contains? (get attrs "rmsk") "repeat")}))

;; ## Linear classifier

(defn- get-vc-inputs
  [attrs normalizer group vc]
  (let [n-vals (normalizer vc)]
    (conj (vec (map #(get n-vals %) attrs)) group)))

(defn- get-train-inputs
  "Retrieve normalized training inputs from VCF file."
  [group in-vcf ctype attrs normalizer ref]
  (let [attr-get (prep-vc-attr-retriever in-vcf ref)]
    (with-open [vcf-iter (get-vcf-iterator in-vcf ref)]
      (->> (parse-vcf vcf-iter)
           (filter #(= ctype (get-classifier-type % attr-get)))
           (map (partial get-vc-inputs attrs normalizer group))
           doall))))

(defn- train-vcf-classifier
  "Do the work of training a variant classifier."
  [ctype attrs base-vcf true-vcf false-vcf ref config]
  (let [normalizer (partial get-vc-attrs-normalized attrs base-vcf ref config)
        inputs (concat (get-train-inputs :a true-vcf ctype attrs
                                         (normalizer true-vcf)
                                         ref)
                       (get-train-inputs :b false-vcf ctype attrs
                                         (normalizer false-vcf)
                                         ref))]
    (when (seq inputs)
      (->> (make-dataset "ds" (conj (vec attrs) {:c [:a :b]}) inputs {:class :c})
           (classifier-train (make-classifier :support-vector-machine :smo))))))

(defn- build-vcf-classifiers
  "Provide a variant classifier based on provided attributes and true/false examples."
  [attrs base-vcf true-vcf false-vcf ref config]
  (letfn [(build-vcf-classifier [ctype]
            (let [out-file (format "%s-%s-classifier.bin" (itx/file-root base-vcf) (ctype-to-str ctype))]
              (if-not (itx/needs-run? out-file)
                (deserialize-from-file out-file)
                (when-let [classifier (train-vcf-classifier ctype attrs base-vcf true-vcf false-vcf
                                                            ref config)]
                  (serialize-to-file classifier out-file)
                  classifier))))]
    (let [ctypes (classifier-types)]
      (zipmap ctypes (map build-vcf-classifier ctypes)))))

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

(defn- vc-passes-w-meta?
  "Check if a variant passes, including external metadata annotations.
   - trusted: pass variants that overlap in the trusted set
   - xspecific: exclude variants specific to a technology or caller
   - otherwise check the variant filtration score, passing those with high scores"
  [vc score meta-getters config]
  (letfn [(meta-has-variants? [kw]
            (has-variants? (get meta-getters kw)
                           (:chr vc) (:start vc) (:end vc)
                           (:ref-allele vc) (:alt-alleles vc)))]
    (cond
     (meta-has-variants? :trusted) true
     (meta-has-variants? :xspecific) false
     :else (> score (get config :min-cscore 0.5)))))

(defn- filter-vc
  "Update a variant context with filter information from classifier."
  [classifiers normalizer attr-get meta-getters config vc]
  (let [attrs (vec (:classifiers config))
        c (get classifiers (get-classifier-type vc attr-get))
        score (-> (make-dataset "ds" (conj attrs :c)
                                 [(get-vc-inputs attrs normalizer -1 vc)]
                                 {:class :c})
                  (make-instance (assoc (normalizer vc) :c -1))
                  (#(classifier-classify c %)))]
    (-> (VariantContextBuilder. (:vc vc))
        (.attributes (assoc (:attributes vc) "CSCORE" score))
        (.filters (when-not (vc-passes-w-meta? vc score meta-getters config)
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
  [base-vcf orig-true-vcf orig-false-vcf meta-files ref config & {:keys [train-w-base?]}]
  (let [out-file (itx/add-file-part base-vcf "cfilter")
        true-vcf (if train-w-base?
                   (get-target-variants base-vcf orig-true-vcf ref "tps")
                   orig-true-vcf)
        false-vcf (if train-w-base?
                    (get-target-variants base-vcf orig-false-vcf ref "fps")
                    orig-false-vcf)
        cs (build-vcf-classifiers (:classifiers config) base-vcf
                                  true-vcf false-vcf ref config)
        normalizer (get-vc-attrs-normalized (:classifiers config) base-vcf ref config base-vcf)
        attr-get (prep-vc-attr-retriever base-vcf ref)]
    (when (itx/needs-run? out-file)
      (println "Filter VCF with" (str cs))
      (with-open [vcf-iter (get-vcf-iterator base-vcf ref)
                  trusted-get (get-vcf-retriever ref (:trusted meta-files))
                  xspecific-get (get-vcf-retriever ref (:xspecific meta-files))]
        (write-vcf-w-template base-vcf {:out out-file}
                              (map (partial filter-vc cs normalizer attr-get
                                            {:trusted trusted-get :xspecific xspecific-get}
                                            config)
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
                             {:trusted (get-train-vcf "trusted")
                              :xspecific (get-train-vcf "xspecific")}
                             (:ref exp) params
                             :train-w-base? (get call :recall false))))
