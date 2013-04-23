(ns bcbio.variation.annotation
  "Annotate variant calls with metrics for assessing false positives
  http://www.broadinstitute.org/gsa/wiki/index.php/VariantAnnotator"
  (:use [bcbio.variation.utils.cgmetrics :only [add-cgmetrics]])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]))

(def ^{:doc "Standard annotations applied to variants"} std-annotations
  ["AlleleBalance" "BaseQualityRankSumTest" "DepthOfCoverage"
   "FisherStrand" "GCContent" "HaplotypeScore" "HomopolymerRun"
   "MappingQualityRankSumTest" "MappingQualityZero"
   "MeanNeighboringBaseQuality" "QualByDepth"
   "ReadPosRankSumTest" "RMSMappingQuality"
   "DepthPerAlleleBySample" "AlleleBalanceConfidenceInterval"
   "MostProbableGenotype" "ReadMeanLen" "ReadMeanPos"
   "ReadPosEndDist" "MinFreeEnergy" "ShannonEntropy"])

(defn add-gatk-annotations
  "Add GATK annotation metrics to variant calls."
  [in-vcf align-bam ref & {:keys [out-dir intervals annos]}]
  {:pre [(not (nil? align-bam))]}
  (let [file-info {:out-vcf (itx/add-file-part in-vcf "annotated" out-dir)}
        ready-annos (if annos annos std-annotations)
        args (concat ["-R" ref
                      "-I" align-bam
                      "--variant" in-vcf
                      "--allow_potentially_misencoded_quality_scores"
                      "-o" :out-vcf]
                     (reduce #(concat %1 ["-A" %2]) [] ready-annos)
                     (broad/gatk-cl-intersect-intervals intervals ref :vcf in-vcf))]
    (broad/index-bam align-bam)
    (broad/run-gatk "VariantAnnotator" args file-info {:out [:out-vcf]})
    (:out-vcf file-info)))

(defn add-variant-annotations
  "Flexible addition of additions to a variant file.
  Handles GATK annotations and Complete Genomics metrics."
  [vcf-file bam-file ref-file call & {:keys [out-dir intervals]}]
  (let [x (get call :annotate "")
        ann (cond
             (true? x) "gatk"
             (false? x) ""
             :else x)]
    (cond
     (and (= ann "gatk") (not (nil? bam-file)))
     (add-gatk-annotations vcf-file bam-file ref-file :out-dir out-dir :intervals intervals)
     (.contains ann "masterVar")
     (add-cgmetrics vcf-file ann ref-file :out-dir out-dir)
     :else vcf-file)))

