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
               fb-vcf (str (fs/file data-dir "freebayes-calls.vcf"))
               c-neg-vcf (str (fs/file data-dir "sv-indels-gatk.vcf"))
               top-out (itx/add-file-part top-vcf "topsubset")
               dip-vcf (str (fs/file data-dir "phasing-input-diploid.vcf"))
               dip-out (itx/add-file-part dip-vcf "haploid")
               c-out (itx/add-file-part top-vcf "cfilter")
               cbin-out (str (itx/file-root top-vcf) "-classifier.bin")]
           (doseq [x (concat [top-out dip-out c-out cbin-out])]
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
  (with-open [vcf-s (get-vcf-source top-vcf ref)]
    (let [vcf-iter (parse-vcf vcf-s)
          attrs ["AD" "QUAL" "DP"]
          normalizer (get-vc-attrs-normalized attrs top-vcf ref)]
      (first (#'bcbio.variation.filter.classify/get-train-inputs
              1 top-vcf attrs normalizer ref)) => (just [0.0 (roughly 0.621) 1.0 1])
      (-> (get-vc-attr-ranges attrs top-vcf ref) (get "DP")) => [193.5 250.0]
      (-> (get-vc-attr-ranges attrs fb-vcf ref) (get "AD")) => (just [0.0
                                                                      (roughly 0.41239)])
      (get-vc-attrs (first vcf-iter) attrs) => {"AD" 0.0 "QUAL" 5826.09 "DP" 250.0}
      (-> (first vcf-iter) normalizer (get "QUAL")) => (roughly 0.621))))

(facts "Final filtration of variants using classifier"
  (filter-vcf-w-classifier top-vcf top-vcf c-neg-vcf ref
                           {:classifiers ["AD" "QUAL" "DP"] :min-cscore 0.5}) => c-out)
