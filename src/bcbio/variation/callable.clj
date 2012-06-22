(ns bcbio.variation.callable
  "Identify callable bases from a BAM alignment file.
  Help differentiate positions where we can not assess variation"
  (:import [org.broad.tribble.bed BEDCodec]
           [org.broad.tribble.index IndexFactory]
           [org.broad.tribble.source BasicFeatureSource])
  (:use [clojure.java.io]
        [bcbio.align.ref :only [sort-bed-file]])
  (:require [clojure.string :as string]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]))

(defn identify-callable
  "Identify callable bases from the provided alignment file."
  [align-bam ref & {:keys [out-dir intervals]}]
  (let [base-dir (if (or (nil? out-dir)
                         (fs/writeable? (fs/parent align-bam)))
                   (fs/parent align-bam)
                   out-dir)
        base-fname (str (file base-dir (-> align-bam fs/base-name itx/file-root)))
        file-info {:out-bed (format "%s-callable.bed" base-fname)
                   :out-summary (format "%s-callable-summary.txt" base-fname)}
        args (concat ["-R" ref
                      "-I" align-bam
                      "--out" :out-bed
                      "--summary" :out-summary]
                     (broad/gatk-cl-intersect-intervals intervals ref))]
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
  [bed-file ref-file]
  (let [batch-size 500
        work-bed (sort-bed-file bed-file ref-file)
        idx (IndexFactory/createIntervalIndex (file work-bed) (BEDCodec.) batch-size)]
    (BasicFeatureSource. work-bed idx (BEDCodec.))))

(defn get-callable-bed
  "Create BED file of callable regions from the BAM alignment file.
  Pass the callable BED to GATK for subsetting based on callable intervals.
  Since the callable output is 1-based inclusive, this converts to 0-based
  intervals on output."
  [align-bam ref & {:keys [out-dir intervals]}]
  (let [orig-bed-file (identify-callable align-bam ref :out-dir out-dir
                                         :intervals intervals)
        out-file (itx/add-file-part orig-bed-file "intervals")]
    (with-open [source (get-bed-source orig-bed-file ref)
                wtr (writer out-file)]
      (doseq [f (.iterator source)]
        (when (= (.getName f) "CALLABLE")
          (.write wtr (format "%s\t%s\t%s\n" (.getChr f)
                              (- (.getStart f) 2) (.getEnd f))))))
    out-file))

;; ## Multiple callables

(defprotocol CallableChecker
  "Provide callable checker for potentially multiple inputs"
  (has-callers? [this])
  (is-callable? [this space start end]))

(defrecord BamCallable [sources check-fn]
  CallableChecker
  (has-callers? [_]
    (not (empty? sources)))
  (is-callable? [_ space start end]
    (letfn [(source-is-callable? [source space start end]
              (if (<= start end)
                (> (count (features-in-region source space start end)) 0)
                false))]
      (if (empty? sources)
        true
        (check-fn #(source-is-callable? % space start end) sources))))
  java.io.Closeable
  (close [_]
    (doseq [x sources]
      (.close x))))

(defn get-callable-checker
  "Retrieve generalized callabilitu checkers that handles multiple file inputs.
  Checks if a chromosome start end region is callable based on reads in input BAM files."
  [bam-files ref & {:keys [out-dir intervals check-fn]
                 :or {check-fn some}}]
  (let [work-bam-files (remove nil? (if (coll? bam-files) bam-files [bam-files]))
        sources (map #(-> (get-callable-bed % ref :out-dir out-dir :intervals intervals)
                          (get-bed-source ref))
                     work-bam-files)]
    (BamCallable. sources check-fn)))
