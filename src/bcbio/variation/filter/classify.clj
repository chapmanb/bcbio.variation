(ns bcbio.variation.filter.classify
  "Provide classification based filtering for variants."
  (:import [org.broadinstitute.variant.variantcontext VariantContextBuilder]
           [org.broadinstitute.variant.vcf VCFHeader VCFInfoHeaderLine
            VCFHeaderLineType VCFFilterHeaderLine VCFHeaderLineCount])
  (:use [ordered.set :only [ordered-set]]
        [clojure.math.combinatorics :only [cartesian-product]]
        [clj-ml.utils :only [serialize-to-file deserialize-from-file]]
        [clj-ml.data :only [make-dataset dataset-set-class make-instance]]
        [clj-ml.classifiers :only [make-classifier classifier-train
                                   classifier-evaluate classifier-classify]]
        [bcbio.variation.filter.util :only [remove-cur-filters]]
        [bcbio.variation.filter.attr :only [get-vc-attrs-normalized prep-vc-attr-retriever]]
        [bcbio.variation.filter.intervals :only [pipeline-combine-intervals]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template
                                               get-vcf-iterator has-variants?
                                               get-vcf-retriever variants-in-region]])
  (:require [clojure.string :as string]
            [me.raynes.fs :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.variation.filter.trusted :as trusted]
            [bcbio.variation.filter.rules :as rules]
            [bcbio.variation.variantcontext :as gvc]))

;; ## Split variants for classification

(defn- classifier-types
  "Define splitting of classifiers based on variant characteristics."
  [attr-key]
  (let [variant-types [:snp :complex]
        zygosity [:hom :het]]
    (map (fn [[vtype z]] {:variant-type vtype
                          :attr-key attr-key
                          :zygosity z})
         (cartesian-product variant-types zygosity))))

(defn- ctype-to-str
  "Convert a classifier types into a string name for output files."
  [x]
  (str (name (:attr-key x)) "-"
       (name (:variant-type x)) "_"
       (name (:zygosity x))))

(defn- get-classifier-type
  "Map variant types to specialized classifiers."
  [vc attr-key attr-get]
  {:variant-type (case (:type vc)
                   "SNP" :snp
                   :complex)
   :attr-key attr-key
   :zygosity (rules/vc-zygosity vc)})

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
           (filter #(= ctype (get-classifier-type % (:attr-key ctype) attr-get)))
           (map (partial get-vc-inputs attrs normalizer group))
           doall))))

(defn- get-dataset
  [attrs inputs]
  (make-dataset "ds" (conj (vec attrs) {:c [:pass :fail]}) inputs {:class :c}))

(defn- train-vcf-classifier
  "Do the work of training a variant classifier."
  [ctype attrs pre-normalizer true-vcf false-vcf ref config]
  (let [config (merge {:normalize :default} config)
        inputs (concat (get-train-inputs :pass true-vcf ctype attrs
                                         (pre-normalizer true-vcf)
                                         ref)
                       (get-train-inputs :fail false-vcf ctype attrs
                                         (pre-normalizer false-vcf)
                                         ref))
        classifier (case (keyword (get config :classifier-type :svm))
                     :svm (make-classifier :support-vector-machine :smo
                                           {:complexity-constant 100.0})
                     :svm-rbf (make-classifier :support-vector-machine :smo
                                               {:kernel-function {:radial-basis {:gamma 0.01}}
                                                :complexity-constant 100000.0})
                     :random-forest (make-classifier :decision-tree :random-forest
                                                     {:num-trees-in-forest 50
                                                      :num-features-to-consider
                                                      (-> attrs count Math/sqrt Math/ceil int)}))]
    (when (seq inputs)
      (->> (get-dataset attrs inputs)
           (classifier-train classifier)))))

(defn- build-vcf-classifiers
  "Provide a variant classifier based on provided attributes and true/false examples."
  [attr-map pre-normalizer base-vcf true-vcf false-vcf ref config out-dir]
  (letfn [(build-vcf-classifier [ctype attrs]
            (let [out-dir (if (nil? out-dir) (str (fs/parent base-vcf)) out-dir)
                  out-file (format "%s/%s-%s-classifier.bin" out-dir
                                   (fs/name base-vcf) (ctype-to-str ctype))]
              (if-not (itx/needs-run? out-file)
                (deserialize-from-file out-file)
                (when-let [classifier (train-vcf-classifier ctype attrs pre-normalizer true-vcf false-vcf
                                                            ref config)]
                  (serialize-to-file classifier out-file)
                  classifier))))]
    (let [ctypes (mapcat classifier-types (keys attr-map))]
      (zipmap ctypes (map #(build-vcf-classifier % (get attr-map (:attr-key %))) ctypes)))))

(defn- add-cfilter-header
  "Add details on the filtering to the VCF file header."
  [attrs]
  (fn [_ header]
    (let [str-attrs (map (fn [[k v]] (str (name k) ": " (string/join "," v))) attrs)
          desc (str "Classification filters based on true/false positives for: "
                    (string/join "; " str-attrs))
          new #{(VCFInfoHeaderLine. "CFILTERS" VCFHeaderLineCount/UNBOUNDED
                                    VCFHeaderLineType/String desc)
                (VCFFilterHeaderLine. "CScoreFilter" "Based on classifcation CFILTERS")}]
      (VCFHeader. (apply ordered-set (concat (.getMetaDataInInputOrder header) new))
                  (.getGenotypeSamples header)))))

(defn- vc-passes-w-meta?
  "Check if a variant passes, including external metadata annotations.
   - trusted: pass variants that overlap in the trusted set
   - xspecific: exclude variants specific to a technology or caller
   - otherwise check the variant filters that failed, passing those that are clean"
  [vc c-filters meta-getters config]
  (letfn [(meta-has-variants? [kw]
            (has-variants? (get meta-getters kw)
                           (:chr vc) (:start vc) (:end vc)
                           (:ref-allele vc) (:alt-alleles vc)))]
    (cond
     (meta-has-variants? :trusted) true
     (meta-has-variants? :xspecific) false
     (empty? c-filters) true
     :else false)))

(defn- filter-vc
  "Update a variant context with filter information from classifier."
  [cs normalizer attr-get meta-getters config vc]
  (letfn [(check-attrgroup-classifier [[attr-key attrs]]
            (let [c (get cs (get-classifier-type vc attr-key attr-get))
                  val (get-vc-inputs attrs normalizer :fail vc)
                  score (classifier-classify c (-> (get-dataset attrs 1)
                                                   (make-instance val)))]
              (when (pos? score) attr-key)))]
    (let [c-filters (->> (:classifiers config)
                         (map check-attrgroup-classifier)
                         (remove nil?))]
      (-> (VariantContextBuilder. (:vc vc))
          (.attributes (assoc (:attributes vc) "CFILTERS" (if (empty? c-filters)
                                                            "None"
                                                            (string/join "," (map name c-filters)))))
          (.filters (when-not (vc-passes-w-meta? vc c-filters meta-getters config)
                      #{"CScoreFilter"}))
          .make))))

(defmulti get-train-variants
  "Retrieve variants to use for true/false positive training.
   Dispatches based on approach used. For recalling, we can pull directly
   from input files"
  (fn [orig-file train-files call exp config call-type out-dir]
    (let [is-recall (get call :recall false)
          recall-approach (keyword (get-in exp [:params :compare-approach] :consensus))]
      (if is-recall
        (if (= :consensus recall-approach)
          [:recall (keyword call-type) (if (:round train-files) :iterate :final)]
          [:recall :rewrite])
        :default))))

(defmethod get-train-variants [:recall :rewrite]
  ^{:doc "Retrieve variants from original file based on variants in target file."}
  [orig-file target-files _ exp _ ext out-dir]
  (letfn [(get-orig-variants [retriever vc]
            (->> (variants-in-region retriever (:chr vc) (:start vc) (:end vc))
                 (filter #(= (:start %) (:start vc)))
                 (map :vc)))]
    (let [out-file (itx/add-file-part orig-file ext out-dir)
          target-file (get target-files (keyword ext))]
      (when (itx/needs-run? out-file)
        (with-open [vcf-iter (get-vcf-iterator target-file (:ref exp))
                    retriever (get-vcf-retriever (:ref exp) orig-file)]
          (write-vcf-w-template orig-file {:out out-file}
                                (mapcat (partial get-orig-variants retriever)
                                        (parse-vcf vcf-iter))
                                (:ref exp))))
      out-file)))

(defmethod get-train-variants [:recall :fps :iterate]
  ^{:doc "Identify false positive variants directly from recalled consensus calls.
          These contain the `set` key value pair with information about supporting
          calls. We filter variants that have low support from multiple callers, then
          compare based on novelty in dbSNP. Novel and know have different inclusion
          parameters derived from examining real true/false calls in replicate
          experiments. The logic for inclusion is:
          - Variants with support from less than `fp-freq` percentage of callers.
            This defaults to less than 25% of callers used.
          - We exclude low mapping quality reads, which end up being non-representative
            of more general cases since they are poorly represented in true positives.
            This is worth looking at more for generalizing the approach to these regions.
          - Include indels in low entropy regions which are often problematic.
          - Include novel variants not found in dbSNP that have low read support.
          - Include known variants, in dbSNP, depending on type:
             - SNP: include SNPs with high likelihood of being ref"}
  [orig-file _ call exp _ ext out-dir]
  (let [passes-rules? (rules/vc-checker orig-file call exp)]
    (letfn [(is-potential-fp? [vc]
              (or (passes-rules? vc
                                 :yes [:below-call-support :high-map-quality :het-snp :low-confidence])
                  (passes-rules? vc
                                 :yes [:below-call-support :high-map-quality :novel])
                  (passes-rules? vc
                                 :yes [:below-call-support :high-map-quality :low-confidence]
                                 :no [:novel])))]
    (gvc/select-variants orig-file is-potential-fp? ext (:ref exp)
                         :out-dir out-dir))))

(defmethod get-train-variants [:recall :tps :iterate]
  ^{:doc "Identify true positive training variants directly from recalled consensus.
          Use variants found in all input callers, then restrict similarly to false
          positives to maintain representative sets. We restrict by lower depth and
          problematic reference likelihoods. We also include high confidence calls
          with lower supporting calls to keep a wider range: these include SNPs with
          a low likelihood of being reference and known indels."}
  [orig-file _ call exp _ ext out-dir]
  (let [passes-rules? (rules/vc-checker orig-file call exp)]
    (letfn [(is-tp? [vc]
              (passes-rules? vc :yes [:all-callers :flex-low-confidence]))]
      (gvc/select-variants orig-file is-tp? ext (:ref exp)
                           :out-dir out-dir))))

(defmethod get-train-variants [:recall :trusted :iterate]
  ^{:doc "Retrieve set of trusted variants based on input parameters and recalled consensus."}
  [orig-file _ call exp params ext out-dir]
  (let [calls (remove #(= (:name %) (:name call)) (:calls exp))]
    (letfn [(is-trusted? [vc]
              (when-let [trusted (:trusted params)]
                (trusted/is-trusted-variant? vc trusted calls)))]
      (gvc/select-variants orig-file is-trusted? ext (:ref exp)
                           :out-dir out-dir))))

(defmethod get-train-variants [:recall :trusted :final]
  [orig-file train-files call exp params ext out-dir]
  (get-train-variants orig-file (assoc train-files :round 1) call exp params ext out-dir))

(defmethod get-train-variants [:recall :xspecific :final]
  ^{:doc "Retrieve specific variants to exclude, handling variants falling below untrusted thresh.
          The `untrusted` keyword in the configuration parameters specifies the threshold to use."}
  [orig-file train-files call exp params ext out-dir]
  (with-open [xspecific-get (gvc/get-vcf-retriever (:ref exp) (:xspecific train-files))]
    (let [calls (remove #(= (:name %) (:name call)) (:calls exp))]
      (letfn [(xspecific? [vc]
                (has-variants? xspecific-get
                               (:chr vc) (:start vc) (:end vc)
                               (:ref-allele vc) (:alt-alleles vc)))
              (untrusted? [vc]
                (when-let [untrusted (:untrusted params)]
                      (not (trusted/is-trusted-variant? vc untrusted calls))))
              (is-untrusted? [vc]
                (or (untrusted? vc) (xspecific? vc)))]
        (gvc/select-variants orig-file is-untrusted? ext (:ref exp)
                             :out-dir out-dir)))))

(defmethod get-train-variants [:recall :tps :final]
  ^{:doc "Iteratively identify true positive variants: low support variants
          that pass the previous round of filtering."}
  [orig-file train-files call exp params ext out-dir]
  (let [passes-rules? (rules/vc-checker orig-file call exp)
        out-file (itx/add-file-part orig-file "tps" out-dir)]
    (letfn [(low-support-novel? [vc]
              (passes-rules? vc
                             :yes [:below-call-support :het-snp :novel :low-depth]))
            (is-previous-tp? [vc]
              (when-not (low-support-novel? vc)
                (or
                 (passes-rules? vc
                                :yes [:below-call-support :passes-filter]
                                :no [:problem-allele-balance
                                     :low-confidence-novel-het-snp])
                 (passes-rules? vc
                                :yes [:below-call-support :het-snp :good-pl]
                                :no [:problem-allele-balance :novel]))))]
      (when (itx/needs-run? out-file)
        (-> (:prev train-files)
            (gvc/select-variants is-previous-tp? ext (:ref exp)
                                 :out-dir out-dir)
            (remove-cur-filters (:ref exp))
            (fs/rename out-file))))
    out-file))

(defmethod get-train-variants [:recall :fps :final]
  ^{:doc "Iteratively identify false positive variants: low support variants
          that fail the previous round of filtering."}
  [orig-file train-files call exp params ext out-dir]
  (let [passes-rules? (rules/vc-checker orig-file call exp)
        out-file (itx/add-file-part orig-file "fps" out-dir)]
    (letfn [(well-supported-known? [vc]
              (passes-rules? vc
                             :yes [:below-call-support :het-snp]
                             :no [:novel :low-depth]))
            (is-previous-fp? [vc]
              (when-not (well-supported-known? vc)
                (or (passes-rules? vc
                                   :yes [:below-call-support :high-map-quality]
                                   :no [:passes-filter])
                    (passes-rules? vc
                                   :yes [:below-call-support :high-map-quality
                                         :het-snp :low-confidence :novel]))))]
      (when (itx/needs-run? out-file)
        (-> (:prev train-files)
            (gvc/select-variants is-previous-fp? ext (:ref exp)
                                 :out-dir out-dir)
            (remove-cur-filters (:ref exp))
            (fs/rename out-file))))
    out-file))

(defmethod get-train-variants :default
  ^{:doc "By default, return the prepped training file with no changes."}
  [_ train-files _ _ _ ext _]
  (get train-files (keyword ext)))

(defn filter-vcf-w-classifier
  "Filter an input VCF file using a trained classifier on true/false variants."
  [base-vcf train-files call exp config]
  (let [out-dir (when-let [tround (:round train-files)]
                  (str (fs/file (fs/parent base-vcf) "trainround") tround))
        out-file (itx/add-file-part base-vcf "cfilter" out-dir)]
    (when (and out-dir (not (fs/exists? out-dir)))
      (fs/mkdirs out-dir))
    (when (itx/needs-run? out-file)
      (let [true-vcf (get-train-variants base-vcf train-files call exp config
                                         "tps" out-dir)
            false-vcf (get-train-variants base-vcf train-files call exp config
                                          "fps" out-dir)
            trusted-vcf (get-train-variants base-vcf train-files call exp
                                            config "trusted" out-dir)
            xspecific-vcf (get-train-variants base-vcf train-files call exp config
                                              "xspecific" out-dir)
            ref (:ref exp)
            pre-normalizer (get-vc-attrs-normalized (apply concat (vals (:classifiers config)))
                                                    base-vcf ref config)
            cs (build-vcf-classifiers (:classifiers config) pre-normalizer base-vcf
                                      true-vcf false-vcf ref config out-dir)
            config (merge {:normalize :default} config)
            attr-get (prep-vc-attr-retriever base-vcf ref)]
        (println "Filter VCF with" (str cs))
        (with-open [vcf-iter (get-vcf-iterator base-vcf ref)
                    trusted-get (get-vcf-retriever ref trusted-vcf)
                    xspecific-get (get-vcf-retriever ref xspecific-vcf)]
          (write-vcf-w-template base-vcf {:out out-file}
                                (map (partial filter-vc cs (pre-normalizer base-vcf) attr-get
                                              {:trusted trusted-get :xspecific xspecific-get}
                                              config)
                                     (parse-vcf vcf-iter))
                                ref :header-update-fn (add-cfilter-header (:classifiers config))))))
    out-file))

(defn pipeline-classify-filter
  "Fit VCF classification-based filtering into analysis pipeline."
  [in-vcf train-info call exp params config]
  (letfn [(get-train-vcf [type]
            (-> (filter #(= type (:name %)) train-info)
                first
                :file))
          (fix-param-classifiers [params]
            (if (map? (:classifiers params))
              params
              (assoc params :classifiers {:all (:classifiers params)})))
          (flatten-param-classifiers [params]
            (assoc params :classifiers
                   {:all (->> (:classifiers params) vals (apply concat) set vec)}))]
    (pipeline-combine-intervals exp config)
    (let [orig-trains {:tps (get-train-vcf "concordant")
                       :fps (get-train-vcf "discordant")
                       :trusted (get-train-vcf "trusted")
                       :xspecific (get-train-vcf "xspecific")}
          params (fix-param-classifiers params)
          x1 (filter-vcf-w-classifier in-vcf (assoc orig-trains :round 1)
                                      call exp (flatten-param-classifiers params))]
      (filter-vcf-w-classifier in-vcf (assoc orig-trains :prev x1) call exp params))))
