(ns bcbio.variation.evaluate
  "Provide high level summary evaluation of variant results, building off GATK VariantEval."
  (:require [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]))

(defn calc-variant-eval-metrics
  "Compare two variant files with GenotypeConcordance in VariantEval"
  [sample vcf1 vcf2 ref & {:keys [out-base intervals]}]
  (let [file-info {:out-eval (str (itx/file-root (if (nil? out-base) vcf1 out-base)) ".eval")}
        args (concat
              ["-R" ref
               "--out" :out-eval
               "--eval" vcf1
               "--comp" vcf2
               "--sample" sample
               "--doNotUseAllStandardModules"
               "--evalModule" "CompOverlap"
               "--evalModule" "CountVariants"
               "--evalModule" "GenotypeConcordance"
               "--evalModule" "TiTvVariantEvaluator"
               "--evalModule" "ValidationReport"
               "--stratificationModule" "Sample"
               "--stratificationModule" "Filter"]
              (broad/gatk-cl-intersect-intervals intervals))]
    (broad/run-gatk "VariantEval" args file-info {:out [:out-eval]})
    (:out-eval file-info)))

(defn- calc-summary-eval-metrics
  "Run VariantEval providing summary information for a VCF file"
  [vcf ref cmp-files]
  (let [file-info {:out-eval (str (itx/file-root vcf) ".eval")}
        args (concat
              ["-R" ref
               "--out" :out-eval
               "--eval" vcf
               "--evalModule" "ThetaVariantEvaluator"
               "--stratificationModule" "AlleleFrequency"
               "--stratificationModule" "Filter"
               "--stratificationModule" "FunctionalClass"
               "--stratificationModule" "IntervalStratification"]
              (flatten (for [[k v] cmp-files]
                         [(format "--comp:%s" (name k)) v])))]
    (broad/run-gatk "VariantEval" args file-info {:out [:out-eval]})
    {:out-eval file-info}))

(defn summary-eval-metrics
  "Provide high level summary metrics of a single variant file."
  [vcf ref cmp-files])
