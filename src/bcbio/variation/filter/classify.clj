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
            [bcbio.run.itx :as itx]
            [bcbio.variation.filter.trusted :as trusted]
            [bcbio.variation.multiple :as multiple]
            [bcbio.variation.variantcontext :as gvc]))

;; ## Split variants for classification

(defn- classifier-types
  "Define splitting of classifiers based on variant characteristics."
  []
  (let [variant-types [:snp :complex]
        repeats [true false]
        zygosity [:hom :het]]
    (map (fn [[vtype rpt z]] {:variant-type vtype
                              :repetitive rpt
                              :zygosity z})
         (cartesian-product variant-types repeats zygosity))))

(defn- ctype-to-str
  "Convert a classifier types into a string name for output files."
  [x]
  (str (name (:variant-type x)) "_"
       (if (:repetitive x) "rpt" "std") "_"
       (name (:zygosity x))))

(defn- get-classifier-type
  "Map variant types to specialized classifiers."
  [vc attr-get]
  (let [attrs (attr-get ["rmsk"] vc)]
    {:variant-type (case (:type vc)
                     "SNP" :snp
                     :complex)
     :repetitive (contains? (get attrs "rmsk") "repeat")
     :zygosity (if (some #(.startsWith (:type %) "HET") (:genotypes vc)) :het :hom)}))

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
  [attrs pre-normalizer base-vcf true-vcf false-vcf ref config]
  (letfn [(build-vcf-classifier [ctype]
            (let [out-file (format "%s-%s-classifier.bin" (itx/file-root base-vcf) (ctype-to-str ctype))]
              (if-not (itx/needs-run? out-file)
                (deserialize-from-file out-file)
                (when-let [classifier (train-vcf-classifier ctype attrs pre-normalizer true-vcf false-vcf
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
     :else (case score
             0.0 true
             1.0 false
             (throw (Exception. (str "Unexpected score: " score)))))))

(defn- filter-vc
  "Update a variant context with filter information from classifier."
  [classifiers normalizer attr-get meta-getters config vc]
  (let [attrs (vec (:classifiers config))
        c (get classifiers (get-classifier-type vc attr-get))
        val (get-vc-inputs attrs normalizer :fail vc)
        score (classifier-classify c (-> (get-dataset attrs 1)
                                         (make-instance val)))]
    (-> (VariantContextBuilder. (:vc vc))
        (.attributes (assoc (:attributes vc) "CSCORE" score))
        (.filters (when-not (vc-passes-w-meta? vc score meta-getters config)
                    #{"CScoreFilter"}))
        .make)))

(defmulti get-train-variants
  "Retrieve variants to use for true/false positive training.
   Dispatches based on approach used. For recalling, we can pull directly
   from input files"
  (fn [orig-file train-file call exp config call-type]
    (let [is-recall (get call :recall false)
          recall-approach (keyword (get-in exp [:params :compare-approach] :consensus))]
      (if is-recall
        (if (= :consensus recall-approach)
          [:recall (keyword call-type)]
          [:recall :rewrite])
        :default))))

(defmethod get-train-variants [:recall :rewrite]
  ^{:doc "Retrieve variants from original file based on variants in target file."}
  [orig-file target-file _ exp _ ext]
  (letfn [(get-orig-variants [retriever vc]
            (->> (variants-in-region retriever (:chr vc) (:start vc) (:end vc))
                 (filter #(= (:start %) (:start vc)))
                 (map :vc)))]
    (let [out-file (itx/add-file-part orig-file ext)]
      (when (itx/needs-run? out-file)
        (with-open [vcf-iter (get-vcf-iterator target-file (:ref exp))
                    retriever (get-vcf-retriever (:ref exp) orig-file)]
          (write-vcf-w-template orig-file {:out out-file}
                                (mapcat (partial get-orig-variants retriever)
                                        (parse-vcf vcf-iter))
                                (:ref exp))))
      out-file)))

(defmethod get-train-variants [:recall :fps]
  ^{:doc "Identify false positive variants directly from recalled consensus calls.
          These contain the `set` key value pair with information about supporting
          calls. We filter variants that have low support from multiple callers, then
          compare based on novelty in dbSNP. Novel and know have different inclusion
          parameters derived from examining real true/false calls in replicate
          experiments. The logic for inclusion is:
          - A variant must be supported by less than `fp-freq` percentage of callers.
            This defaults to less than 25% of callers used.
          - We exclude low mapping quality reads, which end up being non-representative
            of more general cases since they are poorly represented in true positives.
            This is worth looking at more for generalizing the approach to these regions.
          - Include indels in low entropy regions which are often problematic.
          - Include novel variants not found in dbSNP that have low read support.
          - Include known variants, in dbSNP, depending on type:
             - SNP: include SNPs with high likelihood of being ref
             - indel: include indels near the end of reads with high ref likelihood"}
  [orig-file _ call exp _ ext]
  (let [freq (get call :fp-freq 0.25)
        thresh (Math/ceil (* freq (dec (count (:calls exp)))))
        attr-get (prep-vc-attr-retriever orig-file (:ref exp))]
    (letfn [(below-support-thresh? [vc]
              (-> (multiple/get-vc-set-calls vc (:calls exp))
                  (disj (:name call))
                  count
                  (<= thresh)))
            (get-one-attr [attr vc]
              (-> (attr-get [attr] vc) (get attr)))
            (passes-mapping-quality? [vc]
              (let [mq (get-one-attr "MQ" vc)]
                (or (nil? mq)
                    (> mq 50.0))))
            (low-entropy-indel? [vc]
              (let [entropy (get-one-attr "Entropy" vc)]
                (when-not (nil? entropy)
                  (when (not= "SNP" (:type vc))
                    (< entropy 2.6)))))
            (novel-variant? [vc]
              (contains? #{nil "."} (:id vc)))
            (include-novel? [vc]
              (let [attrs (attr-get ["DP"] vc)]
                (when (not-any? nil? (vals attrs))
                  (< (get attrs "DP") 100.0))))
            (include-known? [vc]
              (let [attrs (attr-get ["DP" "PL" "ReadPosEndDist"] vc)]
                (when (not-any? nil? (vals attrs))
                  (if (= "SNP" (:type vc))
                    (> (get attrs "PL") -20.0)
                    (and (> (get attrs "PL") -20.0)
                         (< (get attrs "ReadPosEndDist") 15.0))))))
            (is-potential-fp? [vc]
              (and (below-support-thresh? vc)
                   (passes-mapping-quality? vc)
                   (or (low-entropy-indel? vc)
                       (if (novel-variant? vc)
                         (include-novel? vc)
                         (include-known? vc)))))]
    (gvc/select-variants orig-file is-potential-fp? ext (:ref exp)))))

(defmethod get-train-variants [:recall :tps]
  ^{:doc "Identify true positive training variants directly from recalled consensus.
          Use variants found in all input callers, then restrict similarly to false
          positives to maintain representative sets. We restrict by lower depth and
          problematic reference likelihoods."}
  [orig-file _ call exp _ ext]
  (let [attr-get (prep-vc-attr-retriever orig-file (:ref exp))]
    (letfn [(include-tp? [vc]
              (let [attrs (attr-get ["DP" "PL" "ReadPosEndDist"] vc)]
                (when (not-any? nil? (vals attrs))
                  (and (< (get attrs "DP") 100.0)
                       (> (get attrs "PL") -20.0)))))
            (is-tp? [vc]
              (when-let [set-val (get-in vc [:attributes "set"])]
                (and (= set-val "Intersection")
                     (include-tp? vc))))]
      (gvc/select-variants orig-file is-tp? ext (:ref exp)))))

(defmethod get-train-variants [:recall :trusted]
  ^{:doc "Retrieve set of trusted variants based on input parameters and recalled consensus."}
  [orig-file _ call exp params ext]
  (let [calls (remove #(= (:name %) (:name call)) (:calls exp))]
    (letfn [(is-trusted? [vc]
              (when-let [trusted (:trusted params)]
                (trusted/is-trusted-variant? vc trusted calls)))]
      (gvc/select-variants orig-file is-trusted? ext (:ref exp)))))

(defmethod get-train-variants :default
  ^{:doc "By default, return the prepped training file with no changes."}
  [_ train-file _ _ _ _]
  train-file)

(defn filter-vcf-w-classifier
  "Filter an input VCF file using a trained classifier on true/false variants."
  [base-vcf orig-true-vcf orig-false-vcf meta-files call exp config]
  (let [out-file (itx/add-file-part base-vcf "cfilter")]
    (when (itx/needs-run? out-file)
      (let [true-vcf (get-train-variants base-vcf orig-true-vcf call exp config "tps")
            false-vcf (get-train-variants base-vcf orig-false-vcf call exp config "fps")
            trusted-vcf (get-train-variants base-vcf (:trusted meta-files) call exp
                                            config "trusted")
            ref (:ref exp)
            pre-normalizer (get-vc-attrs-normalized (:classifiers config) base-vcf ref config)
            cs (build-vcf-classifiers (:classifiers config) pre-normalizer base-vcf
                                      true-vcf false-vcf ref config)
            config (merge {:normalize :default} config)
            attr-get (prep-vc-attr-retriever base-vcf ref)]
        (println "Filter VCF with" (str cs))
        (with-open [vcf-iter (get-vcf-iterator base-vcf ref)
                    trusted-get (get-vcf-retriever ref trusted-vcf)
                    xspecific-get (get-vcf-retriever ref (:xspecific meta-files))]
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
                       :file))]
    (pipeline-combine-intervals exp config)
    (filter-vcf-w-classifier in-vcf (get-train-vcf "concordant")
                             (get-train-vcf "discordant")
                             {:trusted (get-train-vcf "trusted")
                              :xspecific (get-train-vcf "xspecific")}
                              call exp params)))
