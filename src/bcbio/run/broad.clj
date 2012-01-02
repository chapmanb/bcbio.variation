;; High level functions to run software from Broad: GATK, Picard

(ns bcbio.run.broad
  (:import [org.broadinstitute.sting.gatk CommandLineGATK]
           [net.sf.samtools SAMFileReader]
           [net.sf.picard.sam BuildBamIndex])
  (:use [clojure.java.io])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]))

(defn run-gatk [program args file-info map-info]
  "Run a GATK commandline in an idempotent file-safe transaction."
  (if (itx/needs-run? (map #(% file-info) (get map-info :out [])))
    (let [std-args ["-T" program "--phone_home" "NO_ET"]]
      (itx/with-tx-files [tx-file-info file-info (get map-info :out [])]
        (CommandLineGATK/start (CommandLineGATK.)
                               (into-array (itx/subs-kw-files
                                            (concat std-args args)
                                            tx-file-info)))))))

(defn index-bam [in-bam]
  "Generate BAM index, skipping if already present."
  (let [index-file (str in-bam ".bai")]
    (if (itx/needs-run? index-file)
      (BuildBamIndex/createIndex (SAMFileReader. (file in-bam)) (file index-file)))
    index-file))
