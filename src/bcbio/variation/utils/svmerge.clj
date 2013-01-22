(ns bcbio.variation.utils.svmerge
  "Merge structural variants into a set of smaller SNP/indel calls.
   Different software detects larger structural variants, requiring
   a final preparation step combining structural and standard calls,
   removing any smaller calls which overlap large insertions and
   deletions.
   This is currently tuned for fosmid merging and reconstruction but
   has knobs to generalize for diploid merging with appropriate phasing
   of variants."
  (:require [clojure.java.io :as io]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.variation.combine :as combine]
            [bcbio.variation.filter.intervals :as intervals]
            [bcbio.variation.normalize :as normalize]
            [bcbio.variation.structural :as structural]
            [bcbio.variation.variantcontext :as gvc]))

(defn- sv->bed
  "Create a BED file of structural variant regions from input VCF."
  [sv-file ref-file]
  (let [out-file (str (itx/file-root sv-file) "-regions.bed")]
    (with-open [wtr (io/writer out-file)]
      (doseq [vc (structural/parse-vcf-sv sv-file ref-file)]
        (when (contains? #{:DEL :INS} (:sv-type vc))
          (.write wtr (format "%s\t%s\t%s\n" (:chr vc) (dec (:start-ci vc)) (:end-ci vc))))))
    out-file))

(defn into-calls
  "Merge structural variants into calls, updating variants and BED regions to assess."
  [call-file region-file sv-file ref-file]
  (let [out-files {:calls (itx/add-file-part call-file "wsvs")
                   :regions (itx/add-file-part region-file "wsvs")}]
    (when (itx/needs-run? (vals out-files))
      (let [sample (first (intervals/get-sample-names call-file))
            svready-file (normalize/prep-vcf sv-file ref-file sample
                                             :config {:prep-sv-genotype true
                                                      :prep-allele-count 1})
            sv-bed (sv->bed svready-file ref-file)
            call-safesv (-> (intervals/select-by-sample sample call-file nil ref-file
                                                        :ext "safesv"
                                                        :exclude-intervals sv-bed))]
        (-> (combine/combine-variants [call-safesv svready-file] ref-file
                                      :merge-type :full)
            (fs/rename (:calls out-files)))
        (-> (intervals/combine-multiple-intervals region-file [] ref-file
                                                  :combine-rule :union
                                                  :more-beds [sv-bed])
            (fs/rename (:regions out-files)))
        ))
    out-files))