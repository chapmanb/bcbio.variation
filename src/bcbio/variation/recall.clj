(ns bcbio.variation.recall
  "Recall batched sets of variants containing no-call regions.
  Combined variant calls from batches contain regions called in
  some samples but not others. The approach:
    - Split sample into called and no-call variants
    - Re-call the no-call variants using the UnifiedGenotyper
    - Merge previously called and re-called into final set."
  (:import [org.broadinstitute.sting.utils.variantcontext
            Genotype VariantContextBuilder GenotypesContext])
  (:use [bcbio.variation.callable :only [get-callable-checker is-callable?]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template get-vcf-source]])
  (:require [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]))

(defn convert-no-calls
  "Convert no-calls into callable reference and real no-calls."
  [in-vcf align-bam ref & {:keys [out-dir intervals num-alleles]}]
  (letfn [(ref-genotype [g vc]
            (doto (-> vc :vc .getGenotypes GenotypesContext/copy)
              (.replace
               (Genotype/modifyAlleles (:genotype g)
                                       (repeat (if (nil? num-alleles)
                                                 (count (:alleles g))
                                                 num-alleles)
                                               (:ref-allele vc))))))
          (maybe-callable-vc [vc call-source]
            {:pre (= 1 (:num-samples vc))}
            (let [g (-> vc :genotypes first)]
              (if (.isNoCall (-> g :alleles first))
                (if (is-callable? call-source (:chr vc) (:start vc) (:end vc))
                  (-> (VariantContextBuilder. (:vc vc))
                      (.genotypes (ref-genotype g vc))
                      (.make))
                  (-> (VariantContextBuilder. (:vc vc))
                      (.filters #{"NotCallable"})
                      (.make)))
                (:vc vc))))
          (convert-vcs [vcf-source call-source]
            (for [vc (parse-vcf vcf-source)]
              [:out (maybe-callable-vc vc call-source)]))]
    (let [out-file (itx/add-file-part in-vcf "wrefs")]
      (when (itx/needs-run? out-file)
        (with-open [in-vcf-s (get-vcf-source in-vcf ref)
                    call-source (get-callable-checker align-bam ref :out-dir out-dir
                                                      :intervals intervals)]
          (write-vcf-w-template in-vcf {:out out-file}
                                (convert-vcs in-vcf-s call-source) ref)))
      out-file)))