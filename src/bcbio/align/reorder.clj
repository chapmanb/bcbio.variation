(ns bcbio.align.reorder
  "Reorder BAM alignment files to a reference dictionary, potentially swapping naming.
  Handles Human hg19 to GRCh37 naming conversions."
  (:import [net.sf.samtools SAMFileReader SAMFileWriterFactory SAMReadGroupRecord])
  (:use [clojure.java.io]
        [bcbio.align.ref :only [get-seq-dict]]
        [bcbio.run.broad :only [index-bam]])
  (:require [bcbio.run.itx :as itx]))

(defn- updated-bam-header
  "Add updated sequence dictionary and run group information to header."
  [in-bam ref-file call exp]
  (letfn [(update-rgs [rgs]
            (if-not (empty? rgs) rgs
                    [(doto (SAMReadGroupRecord. "1")
                       (.setLibrary (:sample exp))
                       (.setPlatform (get call :platform "illumina"))
                       (.setSample (:sample exp))
                       (.setPlatformUnit (:sample exp)))]))]
    (let [read-groups (update-rgs (-> in-bam .getFileHeader .getReadGroups))]
      (doto (-> in-bam .getFileHeader .clone)
        (.setSequenceDictionary (-> ref-file get-seq-dict))
        (.setReadGroups read-groups)))))

(defn chromsome-remap
  "Retrieve order of chromosomes to fetch and mapping to new index."
  [bam-names ref-names]
  (if (every? #(apply = %) (partition 2 (interleave ref-names bam-names))) []
      (println ref-names bam-names)))

(defn reorder-bam
  "Reorder and remap BAM file to match supplied reference file."
  [bam-file ref-file call exp & {:keys [out-dir]}]
  (let [out-file (itx/add-file-part bam-file "reorder" out-dir)]
    (index-bam bam-file)
    (with-open [in-bam (SAMFileReader. (file bam-file))]
      (let [ref-names (map #(.getSequenceName %) (-> ref-file get-seq-dict .getSequences))
            bam-names (map #(.getSequenceName %) (-> in-bam .getFileHeader .getSequenceDictionary
                                                     .getSequences))
            header (updated-bam-header in-bam ref-file call exp)]
        (println header)
        (println ref-names)
        (println bam-names)))))
