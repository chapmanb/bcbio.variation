(ns bcbio.variation.filter.intervals
  "Combined interval lists from filtered variants prepared via multiple calls.
  Multiple call approaches and technologies result in reduced call regions due
  to coverage. These functions manage creation of reduced BED files."
  (:import [org.broadinstitute.sting.utils.interval IntervalUtils
            IntervalMergingRule IntervalSetRule]
           [org.broadinstitute.sting.utils GenomeLocParser]
           [org.broadinstitute.sting.gatk.datasources.reference ReferenceDataSource])
  (:use [clojure.java.io]
        [bcbio.variation.callable :only [get-callable-bed get-bed-iterator]])
  (:require [bcbio.run.itx :as itx]))

(defn- bed-to-intervals
  [bed-file ref-file loc-parser]
  (with-open [bed-iter (get-bed-iterator bed-file ref-file)]
    (doall (map #(.createGenomeLoc loc-parser %) bed-iter))))

(defn- intersection-of-bed-files
  "Generate list of intervals that intersect in all provided BED files."
  [all-beds ref loc-parser]
  (loop [final []
         intervals (map #(bed-to-intervals % ref loc-parser) all-beds)]
      (if (empty? intervals)
        final
        (recur (IntervalUtils/mergeListsBySetOperator final (first intervals)
                                                      IntervalSetRule/INTERSECTION)
               (rest intervals)))))

(defn combine-multiple-intervals
  "Combine intervals from an initial BED and coverage BAM files."
  [initial-bed align-bams ref & {:keys [out-dir]}]
  (let [all-beds (cons initial-bed (map #(get-callable-bed % ref :out-dir out-dir
                                                           :intervals initial-bed)
                                        align-bams))
        loc-parser (GenomeLocParser. (.getReference (ReferenceDataSource. (file ref))))
        out-file (itx/add-file-part initial-bed "multicombine" out-dir)]
    (when (itx/needs-run? out-file)
      (with-open [wtr (writer out-file)]
        (doseq [x (IntervalUtils/sortAndMergeIntervals
                   loc-parser (intersection-of-bed-files all-beds ref loc-parser)
                   IntervalMergingRule/ALL)]
          (.write wtr (format "%s\t%s\t%s\n" (.getContig x) (dec (.getStart x)) (.getStop x))))))
    out-file))