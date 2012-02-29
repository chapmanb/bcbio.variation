(ns bcbio.variation.annotation
  "Annotate variant calls with metrics for assessing false positives
  http://www.broadinstitute.org/gsa/wiki/index.php/VariantAnnotator"
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]))

(defn add-variant-annotations
  "Add GATK annotation metrics to variant calls."
  [in-vcf align-bam ref]
  {:pre [(not (nil? align-bam))]}
  (let [file-info {:out-vcf (itx/add-file-part in-vcf "annotated")}
        annotations ["BaseQualityRankSumTest" "DepthOfCoverage" "FisherStrand"
                     "GCContent" "HaplotypeScore" "HomopolymerRun"
                     "MappingQualityRankSumTest" "MappingQualityZero"
                     "MeanNeighboringBaseQuality" "QualByDepth"
                     "ReadPosRankSumTest" "RMSMappingQuality"]
        args (concat ["-R" ref
                      "-I" align-bam
                      "--variant" in-vcf
                      "-o" :out-vcf]
                     (reduce #(concat %1 ["-A" %2]) [] annotations))]
    (broad/index-bam align-bam)
    (broad/run-gatk "VariantAnnotator" args file-info {:out [:out-vcf]})
    (:out-vcf file-info)))
