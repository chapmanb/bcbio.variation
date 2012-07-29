(ns bcbio.run.broad
  "High level functions to run software from Broad: GATK, Picard"
  (:import [org.broadinstitute.sting.gatk CommandLineGATK]
           [net.sf.samtools SAMFileReader SAMFileReader$ValidationStringency])
  (:use [clojure.java.io]
        [bcbio.align.ref :only [sort-bed-file create-ref-dict]])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]))

(defn- create-ref-dict-gatk
  "Ensure reference dictionary "
  [args]
  (when-let [ref-file (second (drop-while
                               #(not (contains? #{"-R" "--reference_sequence"} %))
                               args))]
    (create-ref-dict ref-file)))

(defn run-gatk
  "Run a GATK commandline in an idempotent file-safe transaction."
  [program args file-info map-info]
  (when (itx/needs-run? (map #(% file-info) (get map-info :out [])))
    (create-ref-dict-gatk args)
    (let [std-args (concat ["-T" program]
                           (when-not (contains? (set args) "--unsafe")
                             ["--unsafe" "LENIENT_VCF_PROCESSING"]))]
      (itx/with-tx-files [tx-file-info file-info (get map-info :out []) [".idx"]]
        (CommandLineGATK/start (CommandLineGATK.)
                               (into-array (map str (itx/subs-kw-files
                                                     (concat std-args args)
                                                     tx-file-info))))))))

(defn index-bam
  "Generate BAM index, skipping if already present."
  [in-bam]
  (let [index-file (str in-bam ".bai")]
    (when (itx/needs-run? index-file)
      (SAMFileReader/setDefaultValidationStringency SAMFileReader$ValidationStringency/LENIENT)
      (BuildBamIndex/createIndex (SAMFileReader. (file in-bam)) (file index-file)))
    index-file))

(defn gatk-cl-intersect-intervals
  "Supply GATK commandline arguments for interval files, merging via intersection."
  [intervals ref-file & {:keys [vcf]}]
  (cond
   (nil? intervals) (if vcf ["-L" vcf] [])
   (coll? intervals) (concat (flatten (map #(list "-L" %)
                                           (map #(sort-bed-file % ref-file) intervals)))
                             ["--interval_set_rule" "INTERSECTION"])
   :else ["-L" (sort-bed-file intervals ref-file)]))
