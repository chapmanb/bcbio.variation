(ns bcbio.variation.callable
  "Identify callable bases from a BAM alignment file.
  Help differentiate positions where we can not assess variation"
  (:import [org.broad.tribble.bed BEDCodec]
           [org.broad.tribble.index IndexFactory]
           [org.broad.tribble.source BasicFeatureSource])
  (:use [clojure.java.io])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]))

(defn identify-callable
  "Identify callable bases from the provided alignment file."
  [align-bam ref & {:keys [out-dir] :or {out-dir nil}}]
  (let [base-dir (if (nil? out-dir) (fs/parent align-bam) out-dir)
        base-fname (str (file base-dir (-> align-bam fs/base-name itx/file-root)))
        file-info {:out-bed (format "%s-callable.bed" base-fname)
                   :out-summary (format "%s-callable-summary.txt" base-fname)}
        args ["-R" ref
              "-I" align-bam
              "--out" :out-bed
              "--summary" :out-summary]]
    (if-not (fs/exists? base-dir)
      (fs/mkdirs base-dir))
    (broad/index-bam align-bam)
    (broad/run-gatk "CallableLoci" args file-info {:out [:out-bed :out-summary]})
    (:out-bed file-info)))

(defn features-in-region [source space start end]
  (for [f (.query source space start end)]
    {:chr (.getChr f)
     :start (.getStart f)
     :end (.getEnd f)
     :name (.getName f)
     :score (.getScore f)
     :strand (.getStrand f)}))

(defn get-bed-source
  "Provide tribble feature source for a BED formatted file."
  [bed-file]
  (let [batch-size 500
        idx (IndexFactory/createIntervalIndex (file bed-file) (BEDCodec.) batch-size)]
    (BasicFeatureSource. bed-file idx (BEDCodec.))))

(defn callable-checker
  "Provide function to check if a region (chromsome start end) is callable.
  Calculates based on reads in input BAM file."
  [align-bam ref & {:keys [out-dir] :or {out-dir nil}}]
  (if (nil? align-bam) [(fn [& _] true) (java.io.StringReader. "")]
      (let [source (get-bed-source (identify-callable align-bam ref :out-dir out-dir))]
        (letfn [(is-callable? [space start end]
                  (> (count (filter #(= (:name %) "CALLABLE")
                                    (features-in-region source space start end)))
                     0))]
          [is-callable? source]))))

(defn get-callable-bed
  "Create BED file of callable regions from the BAM alignment file.
  Pass the callable BED to GATK for subsetting based on callable intervals."
  [align-bam ref & {:keys [out-dir]}]
  (let [orig-bed-file (identify-callable align-bam ref :out-dir out-dir)
        out-file (itx/add-file-part orig-bed-file "intervals")]
    (with-open [source (get-bed-source orig-bed-file)
                wtr (writer out-file)]
      (doseq [f (.iterator source)]
        (when (= (.getName f) "CALLABLE")
          (.write wtr (format "%s\t%s\t%s\n" (.getChr f)
                              (dec (.getStart f)) (.getEnd f))))))
    out-file))
