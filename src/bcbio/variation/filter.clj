;; Filter variant calls according to supplied criteria.

(ns bcbio.variation.filter
  (:use [clojure.string :only [split]])
  (:require [bcbio.run.broad :as broad]
            [bcbio.run.itx :as itx]))

(defn variant-filter [in-vcf jexl-filters ref]
  "Perform hard variant filtering with supplied JEXL expression criteria."
  (letfn [(jexl-args [x]
            ["--filterName" (str (first (split x #"\s+")) "Filter")
             "--filterExpression" x])]
    (let [file-info {:out-vcf (itx/add-file-part in-vcf "filter")}
          args (concat ["-R" ref
                        "--variant" in-vcf
                        "-o" :out-vcf]
                       (flatten (map jexl-args jexl-filters)))]
      (broad/run-gatk "VariantFiltration" args file-info {:out [:out-vcf]})
      (:out-vcf file-info))))
