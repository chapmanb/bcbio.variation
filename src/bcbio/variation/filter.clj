(ns bcbio.variation.filter
  "Filter variant calls according to supplied criteria."
  (:import [org.broadinstitute.sting.utils.variantcontext
            VariantContextBuilder])
  (:use [clojure.string :only [split]]
        [bcbio.variation.filter.classify :only [pipeline-classify-filter]]
        [bcbio.variation.filter.trusted :only [get-support-vcfs get-trusted-variants]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template
                                               get-vcf-iterator]])
  (:require [bcbio.run.broad :as broad]
            [bcbio.run.itx :as itx]))

(defn jexl-from-config [jexl-filters]
  "Retrieve GATK JEXL commandline expressions from filters."
  (letfn [(jexl-args [x]
            ["--filterName" (str (first (split x #"\s+")) "Filter")
             "--filterExpression" x])]
    (flatten (map jexl-args jexl-filters))))

(defn jexl-filters-from-map
  "Convert a map of metrics names and ranges into JEXL filter expressions"
  [filter-map]
  (letfn [(to-jexl [[metric [min max]]]
            (format "%s < %s && %s > %s" metric min metric max))]
    (map to-jexl filter-map)))

(defn variant-filter
  "Perform hard variant filtering with supplied JEXL expression criteria."
  [in-vcf jexl-filters ref]
  (let [file-info {:out-vcf (itx/add-file-part in-vcf "filter")}
        args (concat ["-R" ref
                      "--variant" in-vcf
                      "-o" :out-vcf
                      "-l" "ERROR"]
                      (jexl-from-config jexl-filters))]
    (broad/run-gatk "VariantFiltration" args file-info {:out [:out-vcf]})
    (:out-vcf file-info)))

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

(defn remove-cur-filters
  "Remove any filter information in the supplied file."
  [in-vcf ref]
  (letfn [(remove-vc-filter [vc]
            [:out (-> (VariantContextBuilder. (:vc vc))
                      (.passFilters)
                      (.make))])]
    (let [out-file (itx/add-file-part in-vcf "nofilter")]
      (when (itx/needs-run? out-file)
        (with-open [vcf-iter (get-vcf-iterator in-vcf ref)]
          (write-vcf-w-template in-vcf {:out out-file}
                                (map remove-vc-filter (parse-vcf vcf-iter))
                                ref)))
      out-file)))

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
                    trusted-info {:name "trusted"
                                  :file (when-let [trusted (:trusted params)]
                                          (get-trusted-variants cmps-by-name support trusted
                                                                exp config))}]
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
                                      (pipeline-classify-filter % (cons trusted-info train-info)
                                                                exp params)))))
                    (#(assoc-in % [fkey :name] (format "%s-%s" (get-in % [fkey :name]) "recal")))
                    (assoc-in [fkey :mod] "recal")
                    (assoc :re-compare true))))
            init-target (map vector all-params [:c1 :c2]))))
