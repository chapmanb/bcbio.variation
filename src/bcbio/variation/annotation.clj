;; Annotate variant calls with metrics for assessing false positives
;; http://www.broadinstitute.org/gsa/wiki/index.php/VariantAnnotator

(ns bcbio.variation.annotation
  (:use [bcbio.variation.compare :only [run-gatk add-file-part
                                        needs-run? index-bam]]))

(defn add-variant-annotations [in-vcf align-bam ref]
  "Add GATK annotation metrics to variant calls."
  (let [out-vcf (add-file-part in-vcf "annotated")
        annotations ["BaseQualityRankSumTest" "DepthOfCoverage" "FisherStrand"
                     "GCContent" "HaplotypeScore" "HomopolymerRun"
                     "MappingQualityRankSumTest" "MappingQualityZero"
                     "QualByDepth" "ReadPosRankSumTest" "RMSMappingQuality"]
        args (concat ["-R" ref
                      "-I" align-bam
                      "--variant" in-vcf
                      "-o" out-vcf]
                     (reduce #(concat %1 ["-A" %2]) [] annotations))]
    (if (needs-run? out-vcf)
      (do
        (index-bam align-bam)
        (run-gatk "VariantAnnotator" args)))
    out-vcf))
