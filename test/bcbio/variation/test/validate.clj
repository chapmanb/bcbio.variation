(ns bcbio.variation.test.validate
  "Final variant filtration and preparation for validation."
  (:use [midje.sweet]
        [bcbio.variation.haploid]
        [bcbio.variation.filter.classify]
        [bcbio.variation.validate]
        [bcbio.variation.variantcontext])
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

(facts "Generalized attribute retrieval from variant contexts"
  (-> (get-vc-attr-ranges ["AD" "QUAL" "DP"] top-vcf ref)
      (get "DP")) => [241.5 250.0]
  (with-open [vcf-s (get-vcf-source top-vcf ref)]
    (let [vcf-iter (parse-vcf vcf-s)
          attrs ["AD" "QUAL" "DP"]
          normalizer (get-vc-attrs-normalized attrs top-vcf ref)]
      (get-vc-attrs (first vcf-iter) attrs) => {"AD" 0.0 "QUAL" 5826.09 "DP" 250.0}
      (-> (first vcf-iter) normalizer (get "QUAL")) => (roughly 0.2943))))
