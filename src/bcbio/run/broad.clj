(ns bcbio.run.broad
  "High level functions to run software from Broad: GATK, Picard"
  (:import [org.broadinstitute.sting.gatk CommandLineGATK]
           [net.sf.samtools SAMFileReader SAMFileReader$ValidationStringency]
           [net.sf.picard.sam BuildBamIndex])
  (:use [clojure.java.io])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]))

(defn run-gatk
  "Run a GATK commandline in an idempotent file-safe transaction."
  [program args file-info map-info]
  (if (itx/needs-run? (map #(% file-info) (get map-info :out [])))
    (let [std-args ["-T" program]]
      (itx/with-tx-files [tx-file-info file-info (get map-info :out [])]
        (CommandLineGATK/start (CommandLineGATK.)
                               (into-array (map str (itx/subs-kw-files
                                                     (concat std-args args)
                                                     tx-file-info))))))))

(defn index-bam
  "Generate BAM index, skipping if already present."
  [in-bam]
  (let [index-file (str in-bam ".bai")]
    (if (itx/needs-run? index-file)
      (do
        (SAMFileReader/setDefaultValidationStringency SAMFileReader$ValidationStringency/LENIENT)
        (BuildBamIndex/createIndex (SAMFileReader. (file in-bam)) (file index-file))))
    index-file))
