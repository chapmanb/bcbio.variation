;; Identify callable bases from a BAM alignment file
;; Help differentiate positions where we can not assess variation

(ns bcbio.variation.callable
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]))

(defn identify-callable [align-bam ref]
  "Identify callable bases from the provided alignment file."
  (let [out-file (format "%s-callable.bed" (itx/file-root align-bam))
        summary-file (format "%s-callable-summary.txt" (itx/file-root align-bam))
        args ["-R" ref
              "-I" align-bam
              "--out" out-file
              "--summary" summary-file]]
    (if-not (fs/exists? out-file)
      (do
        (broad/index-bam align-bam)
        (broad/run-gatk "CallableLoci" args)))
    out-file))
