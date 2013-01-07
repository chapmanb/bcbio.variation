(ns bcbio.variation.filter
  "Filter variant calls according to supplied criteria."
  (:use [clojure.string :only [split]]
        [bcbio.variation.filter.attr :only [get-vc-attr prep-vc-attr-retriever]]
        [bcbio.variation.filter.classify :only [pipeline-classify-filter]]
        [bcbio.variation.filter.specific :only [get-x-specific-variants]]
        [bcbio.variation.filter.trusted :only [get-support-vcfs get-trusted-variants]]
        [bcbio.variation.filter.util :only [remove-cur-filters]]
        [bcbio.variation.metrics :only [to-float passes-filter?]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template
                                               get-vcf-iterator write-vcf-from-filter
                                               select-variants]])
  (:require [clojure.set :as set]
            [clojure.string :as string]
            [bcbio.run.broad :as broad]
            [bcbio.run.itx :as itx]))

(defn jexl-from-config [jexl-filters]
  "Retrieve GATK JEXL commandline expressions from filters."
  (letfn [(jexl-args [x]
            ["--filterName" (str (first (split x #"\s+")) "Filter")
             "--filterExpression" x])]
    (flatten (map jexl-args jexl-filters))))

(defn variant-filter
  "Perform hard variant filtering with supplied JEXL expression criteria."
  [in-vcf jexl-filters ref]
  (let [file-info {:out-vcf (itx/add-file-part in-vcf "filter")}
        args (concat ["-R" ref
                      "--variant" in-vcf
                      "-o" :out-vcf
                      "-l" "ERROR"
                      "--unsafe" "ALL" ;"ALLOW_SEQ_DICT_INCOMPATIBILITY"
                      ]
                      (jexl-from-config jexl-filters))]
    (broad/run-gatk "VariantFiltration" args file-info {:out [:out-vcf]})
    (:out-vcf file-info)))

(defn category-variant-filter
  "Perform hard variant filtration handling both range and category metrics"
  [in-vcf metrics ref & {:keys [remove?]}]
  (let [attr-getter (prep-vc-attr-retriever in-vcf ref)]
    (letfn [(infinity-flag? [x]
              (.contains (str x) "Infinity"))
            (in-range? [[orig-min orig-max] x]
              (let [min (if (infinity-flag? orig-min) (- Integer/MAX_VALUE) orig-min)
                    max (if (infinity-flag? orig-max) Integer/MAX_VALUE orig-max)]
                (and (>= x min) (<= x max))))
            (attr-passes? [got want]
              (cond
               (set? want) (not (empty? (set/intersection got want))) 
               (or (vector? want) (list? want)) (in-range? want got)))
            (passes-metrics? [vc]
              (let [attrs (attr-getter (keys metrics) vc)]
                (and (passes-filter? vc) 
                     (every? (fn [[k v]]
                               (attr-passes? (get attrs k) v)) metrics))))
            (range-to-str [k [min max]]
              (cond
               (infinity-flag? min) (format "%s > %.1f" k max)
               (infinity-flag? max) (format "%s < %.1f" k min)
               :else (format "%s not [%.1f %.1f]" k min max)))
            (metric-to-str [[k v]]
              (cond
               (set? v) (format "%s not [%s]" k (string/join "," v))
               (or (vector? v) (list? v)) (range-to-str k v)))]
      (if remove?
        (select-variants in-vcf passes-metrics? "filter" ref)
        (write-vcf-from-filter in-vcf ref "filter"
                               "ManualRanges" (string/join "; " (map metric-to-str metrics))
                               passes-metrics?)))))

(defn variant-format-filter
  "Perform hard filtering base on JEXL expressions on metrics in the Genotype FORMAT field."
  [in-vcf exps ref]
  (letfn [(format-filter [exp]
            (let [[attr op-str str-val] (string/split exp #" ")
                  val (to-float str-val)
                  op (eval (read-string op-str))]
              (fn [vc]
                (when-let [vc-val (get-vc-attr vc [:format attr] {})]
                  (not (op vc-val val))))))
          (format-filter-multi [exps]
            (let [int-filters (map format-filter exps)]
              (fn [vc]
                (every? true? (map #(% vc) int-filters)))))]
    (write-vcf-from-filter in-vcf ref "ffilter"
                           "FormatRanges" (string/join "; " exps)
                           (format-filter-multi exps))))

(defn- variant-recalibration
  "Perform the variant recalibration step with input training VCF files.
  training-vcfs is a list of `{:file vcf-file :name name-to-use :prior probability}`"
  [in-vcf training-vcfs annotations ref & {:keys [lenient]}]
  (let [base-out (itx/file-root in-vcf)
        file-info {:out-recal (str base-out ".recal")
                   :out-tranch (str base-out ".tranches")
                   :out-r (str base-out "-recalplots.R")}
        args (concat ["-R" ref
                      "-input" in-vcf
                      "-recalFile" :out-recal
                      "-tranchesFile" :out-tranch
                      "-rscriptFile" :out-r
                      "--mode" "BOTH"]
                     (if lenient
                       ["--percentBadVariants" "0.05"
                        "--maxGaussians" "4"]
                       ["--percentBadVariants" "0.03"
                        "--maxGaussians" "10"])
                     (flatten (map (fn [x] ["-an" x]) annotations))
                     (flatten (map (fn [x] [(str "-resource:" (:name x)
                                                 ",known=true"
                                                 ",training=true"
                                                 ",truth=" (:truth x)
                                                 ",prior=" (:prior x)
                                                 ",bad=" (:bad x))
                                            (:file x)])
                                   training-vcfs)))]
    (broad/run-gatk "VariantRecalibrator" args file-info {:out [:out-recal :out-tranch]})
    file-info))

(defn- apply-recalibration
  "Apply variant recalibration to input VCF."
  [in-vcf recal-files ref]
  (let [file-info {:out-vcf (itx/add-file-part in-vcf "recalfilter")}
        args ["-R" ref
              "-input" in-vcf
              "--ts_filter_level" "99.0"
              "--mode" "BOTH"
              "-tranchesFile" (:out-tranch recal-files)
              "-recalFile" (:out-recal recal-files)
              "-o" :out-vcf]]
    (broad/run-gatk "ApplyRecalibration" args file-info {:out [:out-vcf]})
    (:out-vcf file-info)))

(defn variant-recal-filter
  "Perform filtration using variant recalibration based on known variations.
  Training-vcfs is a list of true training sites along with associated
  probability and name."
  [in-vcf training-vcfs annotations ref & {:keys [lenient]}]
  (let [recal-files (variant-recalibration in-vcf training-vcfs annotations ref :lenient lenient)]
    (apply-recalibration in-vcf recal-files ref)))

(defn- get-train-info
  "Retrieve training information for GATK recalibration:
   - No support specified: use the target comparison
   - Support specified and a specific comparison pair
   - Support specified as a single target: use target versus all comparison"
  [cmps-by-name support config]
  (let [support-vcfs (get-support-vcfs cmps-by-name support config)]
      [{:file (:true-positives support-vcfs)
        :name "concordant"
        :truth "true"
        :bad "false"
        :prior 10.0}
       {:file (:false-positives support-vcfs)
        :name "discordant"
        :truth "false"
        :bad "true"
        :prior 10.0}]))

(defn pipeline-recalibration
  "Perform variant recalibration and filtration as part of processing pipeline."
  [cmps-by-name finalizer exp config]
  (let [init-target (get cmps-by-name (:target finalizer)
                         (get cmps-by-name (reverse (:target finalizer))))
        all-params (let [x (:params finalizer)] (if (map? x) [x] x))]
    (reduce (fn [target [params fkey]]
              (let [in-vcf (remove-cur-filters (-> target fkey :file) (:ref exp))
                    support (get params :support (:target finalizer))
                    train-info (get-train-info cmps-by-name support config)
                    trusted-info [{:name "trusted"
                                   :file (when-let [trusted (:trusted params)]
                                           (get-trusted-variants cmps-by-name support trusted
                                                                 exp config))}
                                  {:name "xspecific"
                                   :file (when (:xspecific params)
                                           (get-x-specific-variants cmps-by-name support exp config))}]]
                (-> target
                    (assoc-in [fkey :file]
                              (-> in-vcf
                                  (#(if-let [anns (:annotations params)]
                                      (variant-recal-filter % train-info
                                                            anns (:ref exp)
                                                            :lenient (:lenient params))
                                      %))
                                  (#(if-let [hard-filters (:filters params)]
                                      (variant-filter % hard-filters (:ref exp))
                                      %))
                                  (#(if-not (:classifiers params)
                                      %
                                      (pipeline-classify-filter % (concat trusted-info train-info)
                                                                (get target fkey)
                                                                exp params config)))))
                    (#(assoc-in % [fkey :name] (format "%s-%s" (get-in % [fkey :name]) "recal")))
                    (assoc-in [fkey :mod] "recal")
                    (assoc :re-compare true))))
            init-target (map vector all-params [:c1 :c2]))))
