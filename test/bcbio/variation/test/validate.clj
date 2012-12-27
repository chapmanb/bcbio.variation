(ns bcbio.variation.test.validate
  "Final variant filtration and preparation for validation."
  (:use [midje.sweet]
        [bcbio.variation.haploid :exclude [-main]]
        [bcbio.variation.filter.attr]
        [bcbio.variation.filter.classify]
        [bcbio.variation.filter.intervals]
        [bcbio.variation.filter]
        [bcbio.variation.validate]
        [bcbio.variation.variantcontext :exclude [-main]])
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
               c-out-extras (map #(itx/add-file-part top-vcf %) ["fps" "tps" "trusted"])
               cbin-out (str (itx/file-root top-vcf) "-classifier.bin")
               align-bam (str (fs/file data-dir "aligned-reads.bam"))
               region-bed (str (fs/file data-dir "aligned-reads-regions.bed"))
               exclude-bed (str (fs/file data-dir "aligned-reads-regions-exclude.bed"))
               sorted-out (itx/add-file-part exclude-bed "sorted")
               region-multi-out (itx/add-file-part region-bed "multicombine")
               region-out (fs/glob (fs/file data-dir "aligned-reads-callable*"))
               ffilter-out (itx/add-file-part top-vcf "ffilter")]
           (doseq [x (concat [top-out dip-out c-out cbin-out region-multi-out ffilter-out
                              sorted-out]
                             region-out c-out-extras)]
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
  (with-open [vcf-iter (get-vcf-iterator top-vcf ref)]
    (let [vcf-iter (parse-vcf vcf-iter)
          attrs ["AD" "QUAL" "DP"]
          xtra-attrs (conj attrs "gms_illumina")
          config {:normalize "minmax"}
          normalizer ((get-vc-attrs-normalized attrs top-vcf ref config) top-vcf)
          log-normalizer ((get-vc-attrs-normalized attrs top-vcf ref {:normalize "log"
                                                                      :log-attrs ["DP"]})
                          top-vcf)
          ctype {:variant-type :snp :repetitive false :zygosity :hom}]
      (first (#'bcbio.variation.filter.classify/get-train-inputs
              1 top-vcf ctype xtra-attrs normalizer ref)) => (contains [0.0 (roughly 0.621) 1.0 1]
                                                                      :in-any-order :gaps-ok)
      (-> (get-vc-attr-ranges attrs top-vcf ref {}) (get "DP")) => [193.5 250.0]
      (-> (get-vc-attr-ranges attrs fb-vcf ref {}) (get "AD")) => (just [0.0
                                                                      (roughly 0.41239)])
      (get-vc-attrs (first vcf-iter) xtra-attrs {}) => {"gms_illumina" nil
                                                        "AD" 0.0 "QUAL" 5826.09 "DP" 250.0}
      (get-vc-attr (first vcf-iter) [:format "DP"] {}) => 5
      (-> (first vcf-iter) normalizer (get "QUAL")) => (roughly 0.621)
      (-> (first vcf-iter) log-normalizer (get "DP")) => (roughly 5.521))))

(facts "Filter based on genotype FORMAT queries"
  (variant-format-filter top-vcf ["DP < 10"] ref) => ffilter-out)

(facts "Final filtration of variants using classifier"
  (filter-vcf-w-classifier top-vcf top-vcf c-neg-vcf {} {:recall true} {:ref ref}
                           {:classifiers ["AD" "QUAL" "DP" "PL"]
                            :trusted {:total 1.1}}) => c-out)

(facts "Prepare combined interval lists based on filtering criteria"
  (combine-multiple-intervals region-bed [align-bam] ref
                              :exclude-intervals exclude-bed) => region-multi-out)
