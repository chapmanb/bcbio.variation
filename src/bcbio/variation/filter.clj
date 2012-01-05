;; Filter variant calls according to supplied criteria.

(ns bcbio.variation.filter
  (:use [clojure.string :only [split]])
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
