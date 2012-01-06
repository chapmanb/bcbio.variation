;; Filter variant calls according to supplied criteria.

(ns bcbio.variation.filter
  (:import [org.broadinstitute.sting.utils.variantcontext
            VariantContextBuilder])
  (:use [clojure.string :only [split]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template]])
  (:require [bcbio.run.broad :as broad]
            [bcbio.run.itx :as itx]))

(defn variant-filter [in-vcf jexl-filters ref]
  "Perform hard variant filtering with supplied JEXL expression criteria."
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

(defn- variant-recalibration [in-vcf training-vcfs ref]
  "Perform the variant recalibration step with input training VCF files.
  training-vcfs is a list of {:file vcf-file :name name-to-use :prior probability}"
  (let [annotations ["FS" "HaplotypeScore" "HRun" "MQ" "MQRankSum" "ReadPosRankSum" "QD"]
        base-out (itx/file-root in-vcf)
        file-info {:out-recal (str base-out ".recal")
                   :out-tranch (str base-out ".tranches")
                   :out-r (str base-out "-plots.R")}
        args (concat ["-R" ref
                      "-input" in-vcf
                      "-recalFile" :out-recal
                      "-tranchesFile" :out-tranch
                      "-rscriptFile" :out-r]
                     (flatten (map (fn [x] ["-an" x]) annotations))
                     (flatten (map (fn [x] [(str "-resource:" (:name x)
                                                 ",known=true,training=true,truth=true,prior="
                                                 (:prior x))
                                            (:file x)])
                                   training-vcfs)))]
    (broad/run-gatk "VariantRecalibrator" args file-info {:out [:out-recal :out-tranch]})
    file-info))

(defn- apply-recalibration [in-vcf recal-files ref]
  "Apply variant recalibration to input VCF."
  (let [file-info {:out-vcf (itx/add-file-part in-vcf "recalfilter")}
        args ["-R" ref
              "-input" in-vcf
              "--ts_filter_level" "99.0"
              "-tranchesFile" (:out-tranch recal-files)
              "-recalFile" (:out-recal recal-files)
              "-o" :out-vcf]]
    (broad/run-gatk "ApplyRecalibration" args file-info {:out [:out-vcf]})
    (:out-vcf file-info)))

(defn variant-recalibration-filter [in-vcf training-vcfs ref]
  "Perform filtration using variant recalibration based on known variations.
  Training-vcfs is a list of true training sites along with associated
  probability and name."
  (let [recal-files (variant-recalibration in-vcf training-vcfs ref)]
    (apply-recalibration in-vcf recal-files ref)))

(defn remove-cur-filters [in-vcf ref]
  "Remove any filter information in the supplied file."
  (letfn [(remove-vc-filter [vc]
            [:out (-> (VariantContextBuilder. (:vc vc))
                      (.passFilters)
                      (.make))])]
    (let [out-file (itx/add-file-part in-vcf "nofilter")]
      (write-vcf-w-template in-vcf {:out out-file}
                            (map remove-vc-filter (parse-vcf in-vcf))
                            ref)
      out-file)))

(defn pipeline-recalibration [target ref]
  "Perform variant recalibration and filtration as part of processing pipeline."
  (let [in-vcf (remove-cur-filters (-> target :c1 :file) ref)
        train-info [{:file (-> target :c-files first)
                      :name "concordant"
                     :prior 10.0}]]
    (-> target
        (assoc-in [:c1 :file] (variant-recalibration-filter in-vcf train-info ref))
        (assoc-in [:c1 :mod] "recal"))))
