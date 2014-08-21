(ns bcbio.align.reorder
  "Reorder BAM alignment files to a reference dictionary, potentially swapping naming.
  Handles Human hg19 to GRCh37 naming conversions."
  (:import [htsjdk.samtools SAMFileReader SAMFileWriterFactory SAMReadGroupRecord
            SAMTag SAMFileReader$ValidationStringency])
  (:use [clojure.java.io]
        [bcbio.align.ref :only [get-seq-dict]]
        [bcbio.run.broad :only [index-bam]]
        [bcbio.variation.normalize :only [prep-rename-map]])
  (:require [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]))

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

(defn get-new-chr-order
  "Retrieve order of chromosomes to fetch and mapping to new index."
  [bam-names ref-names ref-file]
  (letfn [(get-bam-name-map [bam-names orig-ref-names]
            (let [ref-names (set orig-ref-names)
                  name-remap (prep-rename-map :GRCh37 ref-file)]
              (reduce (fn [coll x]
                        (assoc coll (cond
                                     (contains? ref-names x) x
                                     (contains? name-remap x) (get name-remap x)
                                     :else (throw (Exception. (str "Could not map " x))))
                               x))
                      {} bam-names)))
          (get-index-map [name-map]
            (let [bam-name-map (reduce (fn [coll [x y]] (assoc coll y x))
                                       {} name-map)]
              (reduce (fn [coll [i x]]
                        (assoc coll i (.indexOf ref-names (get bam-name-map x))))
                      {} (map-indexed vector bam-names))))]
    (when-not (every? #(apply = %) (partition 2 (interleave ref-names bam-names)))
      (let [name-map (get-bam-name-map bam-names ref-names)]
        {:names (remove nil? (map #(get name-map %) ref-names))
         :indexes (get-index-map name-map)}))))

(defn bam-read-seq
  "Lazy sequence for BAM reads from a Picard iterator."
  [iter]
  (lazy-seq
   (when (.hasNext iter)
     (cons (.next iter) (bam-read-seq iter)))))

(defn- write-reorder-bam
  "Write reordered BAM file in specified chromosome order."
  [in-bam out-bam chr-order header]
  (let [default-rg-id (-> header .getReadGroups first .getId)]
    (letfn [(update-read [read]
              (let [new-rg-id (if-let [x (.getAttribute read (.name SAMTag/RG))] x
                                      default-rg-id)]
                (doto read
                  (.setHeader header)
                  (.setReferenceIndex (get (:indexes chr-order)
                                           (.getReferenceIndex read) -1))
                  (.setMateReferenceIndex (get (:indexes chr-order)
                                               (.getMateReferenceIndex read) -1))
                  (.setAttribute (.name SAMTag/RG) new-rg-id))))]
      (doseq [cur-chr (:names chr-order)]
        (with-open [iter (.query in-bam cur-chr 0 0 false)]
          (doseq [read (bam-read-seq iter)]
            (.addAlignment out-bam (update-read read)))))
      (with-open [iter (.queryUnmapped in-bam)]
        (doseq [read (bam-read-seq iter)]
          (.addAlignment out-bam (update-read read)))))))

(defn reorder-bam
  "Reorder and remap BAM file to match supplied reference file."
  [bam-file ref-file call exp & {:keys [out-dir]}]
  (let [out-file (fsp/add-file-part bam-file "reorder" out-dir)]
    (when (itx/needs-run? out-file)
      (index-bam bam-file)
      (SAMFileReader/setDefaultValidationStringency SAMFileReader$ValidationStringency/LENIENT)
      (with-open [in-bam (SAMFileReader. (file bam-file))]
        (let [ref-names (map #(.getSequenceName %) (-> ref-file get-seq-dict .getSequences))
              bam-names (map #(.getSequenceName %) (-> in-bam .getFileHeader .getSequenceDictionary
                                                       .getSequences))
              header (updated-bam-header in-bam ref-file call exp)]
          (if-let [chr-order (get-new-chr-order bam-names ref-names ref-file)]
            (do
              (with-open [out-bam (.makeSAMOrBAMWriter (SAMFileWriterFactory.)
                                                       header true (file out-file))]
                (write-reorder-bam in-bam out-bam chr-order header))
              out-file)
            bam-file))))))

(defn -main [bam-file ref-file sample-name]
  (reorder-bam bam-file ref-file {} {:sample sample-name}))
