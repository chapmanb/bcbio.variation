;; Identify callable bases from a BAM alignment file
;; Help differentiate positions where we can not assess variation

(ns bcbio.variation.callable
  (:import [org.broad.tribble.bed BEDCodec]
           [org.broad.tribble.index IndexFactory]
           [org.broad.tribble.source BasicFeatureSource])
  (:use [clojure.java.io])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]))

(defn identify-callable [align-bam ref & {:keys [out-dir] :or {out-dir nil}}]
  "Identify callable bases from the provided alignment file."
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

(defn bed-feature-source [bed-file]
  "Provide tribble feature source for a BED formatted file."
  (let [batch-size 500
        idx (IndexFactory/createIntervalIndex (file bed-file) (BEDCodec.) batch-size)]
    (BasicFeatureSource. bed-file idx (BEDCodec.))))

(defn callable-interval-tree [align-bam ref & {:keys [out-dir] :or {out-dir nil}}]
  "Retrieve an IntervalTree to retrieve information on callability in a region."
  (bed-feature-source (identify-callable align-bam ref :out-dir out-dir)))

(defn callable-checker [align-bam ref & {:keys [out-dir] :or {out-dir nil}}]
  (if (nil? align-bam) (fn [& _] true)
      (let [source (callable-interval-tree align-bam ref :out-dir out-dir)]
        (letfn [(is-callable? [space start end]
                  (> (count (filter #(= (:name %) "CALLABLE")
                                    (features-in-region source space start end)))
                     0))]
          is-callable?))))
