;; High level functions to run software from Broad: GATK, Picard

(ns bcbio.run.broad
  (:import [org.broadinstitute.sting.gatk CommandLineGATK]
           [net.sf.samtools SAMFileReader]
           [net.sf.picard.sam BuildBamIndex])
  (:use [clojure.java.io])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]))

(defn run-gatk [program args]
  (let [std-args ["-T" program "--phone_home" "NO_ET"]]
    (CommandLineGATK/start (CommandLineGATK.)
                           (into-array (concat std-args args)))))

(defn index-bam [in-bam]
  (let [index-file (str in-bam ".bai")]
    (if (itx/needs-run? index-file)
      (BuildBamIndex/createIndex (SAMFileReader. (file in-bam)) (file index-file)))
    index-file))
