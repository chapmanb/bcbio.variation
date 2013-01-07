(ns bcbio.variation.filter.util
  "Provide useful utilities dealing with filtering of variants"
  (:import [org.broadinstitute.sting.utils.variantcontext
            VariantContextBuilder])
(:require [bcbio.run.itx :as itx]
          [bcbio.variation.variantcontext :as gvc]))

(defn remove-cur-filters
  "Remove any filter information in the supplied file."
  [in-vcf ref]
  (letfn [(remove-vc-filter [vc]
            [:out (-> (VariantContextBuilder. (:vc vc))
                      (.passFilters)
                      (.make))])]
    (let [out-file (itx/add-file-part in-vcf "nofilter")]
      (when (itx/needs-run? out-file)
        (with-open [vcf-iter (gvc/get-vcf-iterator in-vcf ref)]
          (gvc/write-vcf-w-template in-vcf {:out out-file}
                                    (map remove-vc-filter (gvc/parse-vcf vcf-iter))
                                    ref)))
      out-file)))
