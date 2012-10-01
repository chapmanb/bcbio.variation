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
        [bcbio.variation.filter.attr :only [get-vc-attrs-normalized]]
        [bcbio.variation.filter.intervals :only [pipeline-combine-intervals]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template
                                               get-vcf-iterator has-variants?
                                               get-vcf-retriever variants-in-region]])
  (:require [clojure.string :as string]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]))

;; ## Split variants for classification

(def ^{:private true
       :doc "Available specialized classifier groups."}
  classifier-types [:snp :complex])

(defn- get-classifier-type
  "Map variant types to specialized classifiers."
  [vc]
  (case (:type vc)
    "SNP" :snp
    :complex))

;; ## Linear classifier

(defn- get-vc-inputs
  [attrs normalizer group vc]
  (let [n-vals (normalizer vc)]
    (conj (vec (map #(get n-vals %) attrs)) group)))

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
     (meta-has-variants? :xspecific) false
     (meta-has-variants? :trusted) true
     :else (> score (get config :min-cscore 0.5)))))

(defn- filter-vc
  "Update a variant context with filter information from classifier."
  [classifiers normalizer meta-getters config vc]
  (let [attrs (vec (:classifiers config))
        score (-> (make-dataset "ds" (conj attrs :c)
                                 [(get-vc-inputs attrs normalizer -1 vc)]
                                 {:class :c})
                  (make-instance (assoc (normalizer vc) :c -1))
                  (#(classifier-classify (get classifiers (get-classifier-type vc)) %)))]
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
        normalizer (get-vc-attrs-normalized (:classifiers config) base-vcf ref config base-vcf)]
    (when (itx/needs-run? out-file)
      (println "Filter VCF with" (str cs))
      (with-open [vcf-iter (get-vcf-iterator base-vcf ref)
                  trusted-get (get-vcf-retriever ref (:trusted meta-files))
                  xspecific-get (get-vcf-retriever ref (:xspecific meta-files))]
        (write-vcf-w-template base-vcf {:out out-file}
                              (map (partial filter-vc cs normalizer
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
