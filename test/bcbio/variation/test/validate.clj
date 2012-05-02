(ns bcbio.variation.test.validate
  "Final variant filtration and preparation for validation."
  (:use [midje.sweet]
        [bcbio.variation.haploid]
        [bcbio.variation.validate])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]))

(background
 (around :facts
         (let [data-dir (str (fs/file "." "test" "data"))
               ref (str (fs/file data-dir "GRCh37.fa"))
               top-vcf (str (fs/file data-dir "gatk-calls.vcf"))
               top-out (itx/add-file-part top-vcf "topsubset")
               dip-vcf (str (fs/file data-dir "phasing-input-diploid.vcf"))
               dip-out (itx/add-file-part dip-vcf "haploid")]
           (doseq [x (concat [top-out dip-out])]
             (itx/remove-path x))
           ?form)))

(let [finalizer {:target "gatk"
                 :params {:validate {:approach "top" :count 5
                                     :top-metric [{:name "QD" "mod" 1}]}}}]
  (facts "Prepare top variants sorted by metrics."
    (get-to-validate top-vcf finalizer ref) => top-out))

(facts "Convert diploid calls into haploid reference variants."
  (diploid-calls-to-haploid dip-vcf ref) => dip-out)
