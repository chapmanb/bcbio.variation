(ns bcbio.align.ref
  "Deal with reference sequences for alignment and variant calling."
  (:import [org.broadinstitute.sting.gatk.datasources.reference ReferenceDataSource]
           [net.sf.picard.reference ReferenceSequenceFileFactory]
           [net.sf.picard.sam CreateSequenceDictionary])
  (:use [clojure.java.io]
        [ordered.map :only [ordered-map]])
  (:require [clojure.string :as string]
            [bcbio.run.itx :as itx]))

(defn create-ref-dict
  [ref-file]
  (let [dict-file (str (itx/file-root ref-file) ".dict")]
    (when (itx/needs-run? dict-file)
      (.instanceMain (CreateSequenceDictionary.)
                     (into-array [(str "r=" ref-file) (str "o=" dict-file)])))
    dict-file))

(defn get-seq-dict
  "Retrieve Picard sequence dictionary from FASTA reference file."
  [ref-file]
  (create-ref-dict ref-file)
  (ReferenceDataSource. (file ref-file))
  (-> ref-file
      file
      ReferenceSequenceFileFactory/getReferenceSequenceFile
      .getSequenceDictionary))

(defn get-seq-name-map
  "Retrieve map of sequence names to index positions in the input reference.
   This is useful for sorting by position."
  [ref-file]
  (reduce (fn [coll [i x]] (assoc coll x i))
          (ordered-map)
          (map-indexed vector
                       (map #(.getSequenceName %) (.getSequences (get-seq-dict ref-file))))))

(defn extract-sequence
  "Retrieve sequence in the provided region from input reference file."
  [ref-file contig start end]
  (let [seq-ref (ReferenceSequenceFileFactory/getReferenceSequenceFile (file ref-file))
        seq-dict (-> seq-ref .getSequenceDictionary)]
    (when (and (contains? (set (map #(.getSequenceName %) (.getSequences seq-dict))) contig)
               (<= end (.getSequenceLength (.getSequence seq-dict contig))))
      (-> seq-ref
          (.getSubsequenceAt contig start end)
          .getBases
          (#(map char %))
          (#(apply str %))))))

(defn sort-bed-file
  "Sort a BED file relative to the input reference"
  [bed-file ref-file]
  (letfn [(process-line [line]
            (let [parts (if (> (count (string/split line #"\t")) 1)
                          (string/split line #"\t")
                          (string/split line #" "))]
              (let [[chr start end] (take 3 parts)]
                [[chr (Integer/parseInt start) (Integer/parseInt end)] line])))
          (ref-sort-fn [ref-file]
            (let [contig-map (get-seq-name-map ref-file)]
              (fn [x]
                (let [sort-vals (first x)]
                  (vec (cons (get contig-map (first sort-vals))
                             (rest sort-vals)))))))]
    (let [out-file (itx/add-file-part bed-file "sorted")]
      (when (itx/needs-run? out-file)
        (with-open [rdr (reader bed-file)
                    wtr (writer out-file)]
          (doseq [[_ line] (sort-by (ref-sort-fn ref-file)
                                    (map process-line (line-seq rdr)))]
            (.write wtr (str line "\n")))))
      out-file)))

