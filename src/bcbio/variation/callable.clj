;; Identify callable bases from a BAM alignment file
;; Help differentiate positions where we can not assess variation

(ns bcbio.variation.callable
  (:use [bcbio.variation.compare :only [run-gatk file-root index-bam]])
  (:require [fs.core :as fs]))

(defn identify-callable [align-bam ref]
  "Identify callable bases from the provided alignment file."
  (let [out-file (format "%s-callable.bed" (file-root align-bam))
        summary-file (format "%s-callable-summary.txt" (file-root align-bam))
        args ["-R" ref
              "-I" align-bam
              "--out" out-file
              "--summary" summary-file]]
    (if-not (fs/exists? out-file)
      (do
        (index-bam align-bam)
        (run-gatk "CallableLoci" args)))
    out-file))
