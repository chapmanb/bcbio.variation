(ns bcbio.variation.filter
  "Filter variant calls according to supplied criteria."
  (:import [org.broadinstitute.sting.utils.variantcontext
            VariantContextBuilder])
  (:use [clojure.string :only [split]]
        [bcbio.variation.multiple :only [multiple-overlap-analysis]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template
                                               get-vcf-source]])
  (:require [bcbio.run.broad :as broad]
            [bcbio.run.itx :as itx]))

(defn variant-filter
  "Perform hard variant filtering with supplied JEXL expression criteria."
  [in-vcf jexl-filters ref]
  (letfn [(jexl-args [x]
            ["--filterName" (str (first (split x #"\s+")) "Filter")
             "--filterExpression" x])]
    (let [file-info {:out-vcf (itx/add-file-part in-vcf "filter")}
          args (concat ["-R" ref
                        "--variant" in-vcf
                        "-o" :out-vcf
                        "-l" "ERROR"]
                       (flatten (map jexl-args jexl-filters)))]
      (broad/run-gatk "VariantFiltration" args file-info {:out [:out-vcf]})
      (:out-vcf file-info))))

(defn- variant-recalibration
  "Perform the variant recalibration step with input training VCF files.
  training-vcfs is a list of `{:file vcf-file :name name-to-use :prior probability}`"
  [in-vcf training-vcfs annotations ref]
  (let [base-out (itx/file-root in-vcf)
        file-info {:out-recal (str base-out ".recal")
                   :out-tranch (str base-out ".tranches")
                   :out-r (str base-out "-recalplots.R")}
        args (concat ["-R" ref
                      "-input" in-vcf
                      "-recalFile" :out-recal
                      "-tranchesFile" :out-tranch
                      "-rscriptFile" :out-r
                      "--percentBadVariants" "0.03"
                      "--maxGaussians" "10"
                      "--mode" "BOTH"]
                     (flatten (map (fn [x] ["-an" x]) annotations))
                     (flatten (map (fn [x] [(str "-resource:" (:name x)
                                                 ",known=true,training=true,"
                                                 "truth=" (:truth x)
                                                 ",prior=" (:prior x))
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
  [in-vcf training-vcfs annotations ref]
  (let [recal-files (variant-recalibration in-vcf training-vcfs annotations ref)]
    (apply-recalibration in-vcf recal-files ref)))

(defn remove-cur-filters
  "Remove any filter information in the supplied file."
  [in-vcf ref]
  (letfn [(remove-vc-filter [vc]
            [:out (-> (VariantContextBuilder. (:vc vc))
                      (.passFilters)
                      (.make))])]
    (let [out-file (itx/add-file-part in-vcf "nofilter")]
      (with-open [vcf-source (get-vcf-source in-vcf ref)]
        (write-vcf-w-template in-vcf {:out out-file}
                              (map remove-vc-filter (parse-vcf vcf-source))
                              ref))
      out-file)))

(defn- get-train-info
  "Retrieve training information for GATK recalibration:
   - No support specified: use the target comparison
   - Support specified and a specific comparison pair
   - Support specified as a single target: use target versus all comparison"
  [cmps-by-name finalizer config]
  (let [support (get finalizer :support (:target finalizer))
        support-vcfs (if (coll? support)
                       (take 2 (-> cmps-by-name (get support) :c-files))
                       (let [x (multiple-overlap-analysis cmps-by-name config support)]
                         [(:true-positives x) (:false-positives x)]))]
      [{:file (first support-vcfs)
        :name "concordant"
        :truth "true"
        :prior 10.0}
       {:file (second support-vcfs)
        :name "discordant"
        :truth "false"
        :prior 10.0}]))

(defn pipeline-recalibration
  "Perform variant recalibration and filtration as part of processing pipeline."
  [cmps-by-name finalizer exp config]
  (let [target (get cmps-by-name (:target finalizer))
        in-vcf (remove-cur-filters (-> target :c1 :file) (:ref exp))
        hard-filters (get-in finalizer [:params :filters])
        anns (get-in finalizer [:params :annotations])
        train-info (get-train-info cmps-by-name finalizer config)]
    (-> target
        (assoc-in [:c1 :file] (-> in-vcf
                                  (#(if-not anns % (variant-recal-filter % train-info
                                                                         anns (:ref exp))))
                                  (#(if-not hard-filters % (variant-filter % hard-filters
                                                                           (:ref exp))))))
        (#(assoc-in % [:c1 :name] (format "%s-%s" (get-in % [:c1 :name]) "recal")))
        (assoc-in [:c1 :mod] "recal")
        (assoc :re-compare true))))
