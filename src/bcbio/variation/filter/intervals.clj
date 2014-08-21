(ns bcbio.variation.filter.intervals
  "Combined interval lists from filtered variants prepared via multiple calls.
  Multiple call approaches and technologies result in reduced call regions due
  to coverage. These functions manage creation of reduced BED files."
  (:import [org.broadinstitute.gatk.utils.interval IntervalUtils
            IntervalMergingRule IntervalSetRule]
           [org.broadinstitute.gatk.utils GenomeLocParser
            GenomeLocSortedSet]
           [org.broadinstitute.gatk.utils.exceptions UserException$BadInput])
  (:use [clojure.java.io]
        [clojure.set :only [intersection]]
        [bcbio.align.ref :only [get-seq-dict]]
        [bcbio.variation.callable :only [get-callable-bed get-bed-iterator]]
        [bcbio.variation.variantcontext :only [get-vcf-header]])
  (:require [me.raynes.fs :as fs]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]))

;; ## interval VCF subsetting by BED

(defn get-sample-names
  "Retrieve samples identified in the input VCF file."
  [in-vcf]
  (-> in-vcf get-vcf-header .getGenotypeSamples vec))

(defn vcf-sample-name
  "Retrieve the sample name in a provided VCF file, allowing for partial matches."
  [sample in-vcf ref-file]
  (letfn [(sample-match [x choices]
            (let [do-match (filter #(when (.contains % x) %) choices)]
              (when (= 1 (count do-match))
                (first do-match))))]
    (let [vcf-samples (-> in-vcf get-vcf-header .getGenotypeSamples set)]
      (cond
       (contains? vcf-samples sample) sample
       (= 1 (count vcf-samples)) (first vcf-samples)
       :else (sample-match sample vcf-samples)))))

(defn select-by-sample
  "Select only the sample of interest from input VCF files."
  [sample in-file name ref & {:keys [out-dir intervals remove-refcalls ext
                                     exclude-intervals]
                              :or {remove-refcalls false}}]
  (let [base-dir (if (nil? out-dir) (fs/parent in-file) out-dir)
        file-info {:out-vcf (if ext (fsp/add-file-part in-file ext out-dir)
                                (str (fs/file base-dir
                                              (format "%s-%s.vcf" sample name))))}
        args (concat ["-R" ref
                      "--sample_name" (vcf-sample-name sample in-file ref)
                      "--variant" in-file
                      "--unsafe" "ALL" ; "ALLOW_SEQ_DICT_INCOMPATIBILITY"
                      "--out" :out-vcf]
                     (when remove-refcalls ["--excludeNonVariants" "--excludeFiltered"])
                     (when exclude-intervals ["--excludeIntervals" exclude-intervals])
                     (broad/gatk-cl-intersect-intervals intervals ref))]
    (if-not (fs/exists? base-dir)
      (fs/mkdirs base-dir))
    (broad/run-gatk "SelectVariants" args file-info {:out [:out-vcf]})
    (:out-vcf file-info)))

;; ## BED manipulation

(defn- bed-to-intervals
  [bed-file ref-file loc-parser]
  (with-open [bed-iter (get-bed-iterator bed-file ref-file)]
    (doall (map #(.createGenomeLoc loc-parser %) bed-iter))))

(defn- intersect-by-contig
  "Intersect a group of intervals present on a contig."
  [start-intervals combine-rule]
  (loop [final []
         intervals start-intervals]
    (if (empty? intervals)
      final
      (recur (try (IntervalUtils/mergeListsBySetOperator final (first intervals)
                                                         (if (= :union combine-rule)
                                                           IntervalSetRule/UNION
                                                           IntervalSetRule/INTERSECTION))
                  (catch UserException$BadInput e []))
             (rest intervals)))))

;; ## Merge BED files

(defn- prep-intervals-by-contig
  "Intersect and exclude intervals on a contig."
  [start-intervals exclude-intervals loc-parser combine-rule]
  (let [overlaps (intersect-by-contig start-intervals combine-rule)]
    (if (empty? exclude-intervals)
      overlaps
      (let [clean-intervals (->> (group-by #(.getStart %) exclude-intervals)
                                 vals
                                 (map (fn [xs] (sort-by #(.size %) > xs)))
                                 (map first))]
        (-> (GenomeLocSortedSet/createSetFromList loc-parser overlaps)
            (.subtractRegions (GenomeLocSortedSet/createSetFromList loc-parser clean-intervals))
            .toList)))))

(defn intersection-of-bed-files
  "Generate list of intervals that intersect in all provided BED files."
  [all-beds ref loc-parser & {:keys [exclude-bed combine-rule]}]
  (letfn [(intervals-by-chrom [bed-file]
            (group-by #(.getContig %) (bed-to-intervals bed-file ref loc-parser)))
          (get-by-contig [interval-groups contig]
            (map #(get % contig []) interval-groups))]
    (let [interval-groups (map intervals-by-chrom all-beds)
          exclude-by-contig (if exclude-bed (intervals-by-chrom exclude-bed) {})
          contigs (vec (apply intersection (map #(set (keys %)) interval-groups)))]
      (mapcat #(prep-intervals-by-contig (get-by-contig interval-groups %)
                                         (get exclude-by-contig % []) loc-parser
                                         combine-rule)
              contigs))))

(defn combine-multiple
  "Combine intervals from an initial BED and coverage BAM files."
  [initial-bed align-bams ref & {:keys [out-dir name exclude-intervals combine-rule
                                        more-beds]}]
  (let [all-beds (concat [initial-bed] more-beds
                         (map #(get-callable-bed % ref :out-dir out-dir
                                                 :intervals initial-bed)
                              align-bams))
        loc-parser (GenomeLocParser. (get-seq-dict ref))
        out-file (fsp/add-file-part initial-bed
                                    (str (if name (str name "-") "") "multicombine")
                                    out-dir)]
    (when (itx/needs-run? out-file)
      (with-open [wtr (writer out-file)]
        (doseq [x (IntervalUtils/sortAndMergeIntervals
                   loc-parser (intersection-of-bed-files all-beds ref loc-parser
                                                         :exclude-bed exclude-intervals
                                                         :combine-rule combine-rule)
                   IntervalMergingRule/ALL)]
          (.write wtr (format "%s\t%s\t%s\n" (.getContig x) (dec (.getStart x)) (.getStop x))))))
    out-file))

(defn pipeline-combine-intervals
  "Combine multiple intervals as part of processing and filtering pipeline."
  [exp config]
  (let [base-intervals (:intervals exp)
        all-aligns (set (remove nil? (map :align (cons exp (:calls exp)))))]
    (when (and base-intervals (seq all-aligns))
      (combine-multiple base-intervals all-aligns
                        (:ref exp)
                        :exclude-intervals (:exclude-intervals exp)
                        :name (:sample exp)
                        :out-dir (get-in config [:dir :prep] (get-in config [:dir :out]))))))
