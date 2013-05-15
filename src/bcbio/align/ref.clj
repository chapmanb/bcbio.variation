(ns bcbio.align.ref
  "Deal with reference sequences for alignment and variant calling."
  (:import [org.broadinstitute.sting.gatk.datasources.reference ReferenceDataSource]
           [net.sf.picard.reference ReferenceSequenceFileFactory]
           [net.sf.picard.sam CreateSequenceDictionary])
  (:use [clojure.java.io]
        [ordered.map :only [ordered-map]])
  (:require [clojure.string :as string]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]))

(defn create-ref-dict
  [ref-file]
  (let [dict-file (str (itx/file-root ref-file) ".dict")]
    (when (itx/needs-run? dict-file)
      (.instanceMain (CreateSequenceDictionary.)
                     (into-array [(str "r=" ref-file) (str "o=" dict-file)])))
    dict-file))

(defn get-seq-dict*
  "Retrieve Picard sequence dictionary from FASTA reference file."
  [ref-file]
  (create-ref-dict ref-file)
  (ReferenceDataSource. (file ref-file))
  (-> (file ref-file)
      ReferenceSequenceFileFactory/getReferenceSequenceFile
      .getSequenceDictionary))

(def get-seq-dict (memoize get-seq-dict*))

(defn get-seq-dict-and-ref*
  "Retrieve Picard sequence dictionary and reference from FASTA file."
  [ref-file]
  (create-ref-dict ref-file)
  (ReferenceDataSource. (file ref-file))
  (let [seq-ref (ReferenceSequenceFileFactory/getReferenceSequenceFile (file ref-file))
        seq-dict (-> seq-ref .getSequenceDictionary)]
    [seq-dict seq-ref]))

(def get-seq-dict-and-ref (memoize get-seq-dict-and-ref*))

(defn get-seq-name-map
  "Retrieve map of sequence names to index positions in the input reference.
   This is useful for sorting by position."
  [ref-file]
  (reduce (fn [coll [i x]] (assoc coll x i))
          (ordered-map)
          (map-indexed vector
                       (map #(.getSequenceName %) (.getSequences (get-seq-dict ref-file))))))

(defn extract-sequence
  "Retrieve sequence in the provided region from input reference file.
   start and end are 1-based inclusive coordinates (VCF style)"
  [ref-file contig start end]
  (let [[seq-dict seq-ref] (get-seq-dict-and-ref ref-file)]
    (when (and (contains? (set (map #(.getSequenceName %) (.getSequences seq-dict))) contig)
               (<= end (.getSequenceLength (.getSequence seq-dict contig))))
      (-> seq-ref
          (.getSubsequenceAt contig start end)
          .getBases
          (#(map char %))
          (#(apply str %))))))

(defn sort-bed-file
  "Sort a BED file relative to the input reference.
   Takes a IO intensive approach over memory intensive by sorting in blocks
   of chromosomes. `same-time-chrs` handles the tradeoff between speed and
   memory by determining how many chromosomes to process simultaneously."
  [bed-file ref-file]
  (letfn [(process-line [cur-chrs line]
            (let [tab-parts (string/split line #"\t")
                  parts (if (> (count tab-parts) 1)
                          tab-parts
                          (string/split line #" "))]
              (let [[chr start end] (take 3 parts)]
                (when (or (nil? cur-chrs) (contains? cur-chrs chr))
                  [[chr (Integer/parseInt start) (Integer/parseInt end)] line]))))
          (get-input-chrs [bed-file]
            (with-open [rdr (reader bed-file)]
              (->> (line-seq rdr)
                   (map (partial process-line nil))
                   (map ffirst)
                   set)))]
    (let [out-file (itx/add-file-part bed-file "sorted")
          input-chrs (get-input-chrs bed-file)
          same-time-chrs 5]
      (when (or (itx/needs-run? out-file)
                (> (fs/mod-time bed-file) (fs/mod-time out-file)))
        (itx/with-tx-file [tx-out out-file]
          (with-open [wtr (writer tx-out)]
            (doseq [cur-chrs (->> (get-seq-dict ref-file)
                                  .getSequences
                                  (map #(.getSequenceName %))
                                  (filter input-chrs)
                                  (partition-all same-time-chrs)
                                  (map set))]
              (with-open [rdr (reader bed-file)]
                (doseq [[_ line] (->> (line-seq rdr)
                                      (map (partial process-line cur-chrs))
                                      (remove nil?)
                                      (sort-by first))]
                  (.write wtr (str line "\n"))))))))
      out-file)))

