;; Annotate variant calls with metrics for assessing false positives
;; http://www.broadinstitute.org/gsa/wiki/index.php/VariantAnnotator

(ns bcbio.variation.annotation
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]))

(defn add-variant-annotations [in-vcf align-bam ref]
  "Add GATK annotation metrics to variant calls."
  (let [out-vcf (itx/add-file-part in-vcf "annotated")
        annotations ["BaseQualityRankSumTest" "DepthOfCoverage" "FisherStrand"
                     "GCContent" "HaplotypeScore" "HomopolymerRun"
                     "MappingQualityRankSumTest" "MappingQualityZero"
                     "QualByDepth" "ReadPosRankSumTest" "RMSMappingQuality"]
        args (concat ["-R" ref
                      "-I" align-bam
                      "--variant" in-vcf
                      "-o" out-vcf]
                     (reduce #(concat %1 ["-A" %2]) [] annotations))]
    (if (itx/needs-run? out-vcf)
      (do
        (broad/index-bam align-bam)
        (broad/run-gatk "VariantAnnotator" args)))
    out-vcf))
