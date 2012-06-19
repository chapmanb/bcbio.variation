(ns bcbio.variation.callable
  "Identify callable bases from a BAM alignment file.
  Help differentiate positions where we can not assess variation"
  (:import [org.broad.tribble.bed BEDCodec]
           [org.broad.tribble.index IndexFactory]
           [org.broad.tribble.source BasicFeatureSource])
  (:use [clojure.java.io])
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
                     (broad/gatk-cl-intersect-intervals intervals))]
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

(defn sort-bed-file
  [bed-file]
  (letfn [(process-line [line]
            (let [parts (if (> (count (string/split line #"\t")) 1)
                          (string/split line #"\t")
                          (string/split line #" "))]
              (let [[chr start end] (take 3 parts)]
                [[chr (Integer/parseInt start) (Integer/parseInt end)] line])))]
    (let [out-file (itx/add-file-part bed-file "sorted")]
      (when (itx/needs-run? out-file)
        (with-open [rdr (reader bed-file)
                    wtr (writer out-file)]
          (doseq [[_ line] (sort (map process-line (line-seq rdr)))]
            (.write wtr (str line "\n")))))
      out-file)))

(defn get-bed-source
  "Provide tribble feature source for a BED formatted file."
  [bed-file]
  (let [batch-size 500
        work-bed (sort-bed-file bed-file)
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
    (with-open [source (get-bed-source orig-bed-file)
                wtr (writer out-file)]
      (doseq [f (.iterator source)]
        (when (= (.getName f) "CALLABLE")
          (.write wtr (format "%s\t%s\t%s\n" (.getChr f)
                              (- (.getStart f) 2) (.getEnd f))))))
    out-file))

(defn callable-checker
  "Provide function to check if a region (chromsome start end) is callable.
  Calculates based on reads in input BAM file."
  [align-bam ref & {:keys [out-dir intervals]}]
  (if (nil? align-bam) [(fn [& _] true) (java.io.StringReader. "")]
      (let [source (get-bed-source (get-callable-bed align-bam ref :out-dir out-dir
                                                     :intervals intervals))]
        (letfn [(is-callable? [space start end]
                  (if (<= start end)
                    (> (count (features-in-region source space start end)) 0)
                    false))]
          [is-callable? source]))))

;; ## Multiple callables

(defprotocol CloseableCallable
  "Provide callable checker for potentially multiple inputs"
  (has-callers? [this])
  (is-callable? [this space start end])
  (close [this]))

(defrecord BamCallable [callables sources check-fn]
  CloseableCallable
  (has-callers? [_]
    (not (empty? callables)))
  (is-callable? [_ space start end]
    (check-fn #(% space start end) callables))
  (close [_]
    (doseq [x sources]
      (.close x))))

(defn check-any-callable
  "High level checker if any reads in a set of targets have callable bases."
  [target-cmps ref out-dir]
  (let [other-bams (reduce (fn [coll x] (conj coll x))
                           #{} (->> (vals target-cmps)
                                    (map (juxt :c1 :c2))
                                    flatten
                                    (map :align)
                                    (remove nil?)))
        checkers (map #(callable-checker % ref :out-dir out-dir)
                      other-bams)]
    (BamCallable. (map first checkers) (map second checkers) some)))
