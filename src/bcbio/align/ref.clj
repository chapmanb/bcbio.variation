(ns bcbio.align.ref
  "Deal with reference sequences for alignment and variant calling."
  (:import [org.broadinstitute.sting.gatk.datasources.reference ReferenceDataSource]
           [net.sf.picard.reference ReferenceSequenceFileFactory]
           [net.sf.picard.sam CreateSequenceDictionary])
  (:use [clojure.java.io]
        [ordered.map :only [ordered-map]])
  (:require [clojure.string :as string]
            [clojure.java.shell :as shell]
            [iota]
            [me.raynes.fs :as fs]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]))

(defn create-ref-dict
  "Create reference dictionaries required by GATK and Picard.
   Requires samtools command to create *.fai if missing, since
   code to create these is no longer present in GATK."
  [ref-file]
  (let [dict-file (str (fsp/file-root ref-file) ".dict")
        fai-file (str ref-file ".fai")]
    (when (itx/needs-run? dict-file)
      (.instanceMain (CreateSequenceDictionary.)
                     (into-array [(str "r=" ref-file) (str "o=" dict-file)])))
    (when (itx/needs-run? fai-file)
      (shell/sh "samtools" "faidx" ref-file))
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

(defn- prep-bedline-sort
  "Convert a line in a BED file into sort coordinates"
  [bed-file ref-file]
  (let [ref-map (get-seq-name-map ref-file)
        is-tab? (with-open [rdr (reader bed-file)]
                  (.contains (first (drop-while #(.startsWith % "track") (line-seq rdr))) "\t"))]
    (fn [line]
      (when-not (.startsWith line "track")
        (let [parts (if is-tab?
                      (string/split line #"\t")
                      (string/split line #" "))]
          (let [[chr start end] (take 3 parts)]
            [(get ref-map chr) (Integer/parseInt start) (Integer/parseInt end)]))))))

(defn sort-bed-file
  "Sort a BED file relative to the input reference.
   Uses memory mapped indexed files to avoid high memory requirements."
  [bed-file ref-file]
  (let [out-file (fsp/add-file-part bed-file "sorted")]
    (when (or (itx/needs-run? out-file)
              (> (fs/mod-time bed-file) (fs/mod-time out-file)))
      (let [bedline->sort (prep-bedline-sort bed-file ref-file)
            bed-vec (iota/vec bed-file)]
        (itx/with-tx-file [tx-out out-file]
          (with-open [wtr (writer tx-out)]
            (doseq [idx (sort-by #(bedline->sort (nth bed-vec %))
                                 (range (count bed-vec)))]
              (let [outl (nth bed-vec idx)]
                (when-not (.startsWith outl "track")
                  (.write wtr (str outl "\n")))))))))
    out-file))
