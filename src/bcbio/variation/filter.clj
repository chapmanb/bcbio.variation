
(ns bcbio.variation.filter
  "Filter variant calls according to supplied criteria."
  (:import [org.broadinstitute.sting.utils.variantcontext
            VariantContextBuilder])
  (:use [clojure.string :only [split]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template]])
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
  [in-vcf training-vcfs ref]
  (let [annotations ["MQRankSum" "ReadPosRankSum"]
        base-out (itx/file-root in-vcf)
        file-info {:out-recal (str base-out ".recal")
                   :out-tranch (str base-out ".tranches")
                   :out-r (str base-out "-recalplots.R")}
        args (concat ["-R" ref
                      "-input" in-vcf
                      "-recalFile" :out-recal
                      "-tranchesFile" :out-tranch
                      "-rscriptFile" :out-r
                      "--percentBadVariants" "0.03"
                      "--maxGaussians" "4"]
                     (flatten (map (fn [x] ["-an" x]) annotations))
                     (flatten (map (fn [x] [(str "-resource:" (:name x)
                                                 ",known=true,training=true,truth=true,prior="
                                                 (:prior x))
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
              "-tranchesFile" (:out-tranch recal-files)
              "-recalFile" (:out-recal recal-files)
              "-o" :out-vcf]]
    (broad/run-gatk "ApplyRecalibration" args file-info {:out [:out-vcf]})
    (:out-vcf file-info)))

(defn variant-recalibration-filter
  "Perform filtration using variant recalibration based on known variations.
  Training-vcfs is a list of true training sites along with associated
  probability and name."
  [in-vcf training-vcfs ref]
  (let [recal-files (variant-recalibration in-vcf training-vcfs ref)]
    (apply-recalibration in-vcf recal-files ref)))

(defn remove-cur-filters
  "Remove any filter information in the supplied file."
  [in-vcf ref]
  (letfn [(remove-vc-filter [vc]
            [:out (-> (VariantContextBuilder. (:vc vc))
                      (.passFilters)
                      (.make))])]
    (let [out-file (itx/add-file-part in-vcf "nofilter")]
      (write-vcf-w-template in-vcf {:out out-file}
                            (map remove-vc-filter (parse-vcf in-vcf))
                            ref)
      out-file)))

(defn pipeline-recalibration
  "Perform variant recalibration and filtration as part of processing pipeline."
  [target support ref]
  (let [in-vcf (remove-cur-filters (-> target :c1 :file) ref)
        hard-filters ["FS > 60.0" "HRun > 5.0" "QD < 2.0" "HaplotypeScore > 13.0"]
        train-info [{:file (-> support :c-files first)
                      :name "concordant"
                     :prior 10.0}]]
    (-> target
        (assoc-in [:c1 :file] (-> in-vcf
                                  (variant-recalibration-filter train-info ref)
                                  (variant-filter hard-filters ref)))
        (assoc-in [:c1 :mod] "recal"))))
