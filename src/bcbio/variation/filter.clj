;; Filter variant calls according to supplied criteria.

(ns bcbio.variation.filter
  (:use [clojure.string :only [split]])
  (:require [bcbio.run.broad :as broad]
            [bcbio.run.itx :as itx]))

(defn variant-filter [in-vcf filters ref]
  "Perform hard variant filtering with supplied JEXL expression criteria."
  (letfn [(filter-args [f]
            ["--filterName" (str (first (split f #" " 1)) "Filter")
             "--filterExpression" f])]
    (let [file-info {:out-vcf (itx/add-file-part in-vcf "filter")}
          args (concat ["-R" ref
                        "--variant" in-vcf
                        "-o" :out-vct]
                       (map filter-args filters))]
      (broad/run-gatk "VariantFiltration" args file-info {:out [:out-vcf]}))))
