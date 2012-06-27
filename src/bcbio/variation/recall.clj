(ns bcbio.variation.recall
  "Recall batched sets of variants containing no-call regions.
  Combined variant calls from batches contain regions called in
  some samples but not others. The approach:
    - Split sample into called and no-call variants
    - Re-call the no-call variants using the UnifiedGenotyper
    - Merge previously called and re-called into final set.
  http://www.broadinstitute.org/gsa/wiki/index.php/Merging_batched_call_sets"
  (:import [org.broadinstitute.sting.utils.variantcontext
            Genotype VariantContextBuilder GenotypesContext]
           [org.broadinstitute.sting.utils.codecs.vcf VCFHeader])
  (:use [ordered.set :only [ordered-set]]
        [bcbio.variation.annotation :only [std-annotations]]
        [bcbio.variation.callable :only [get-callable-checker is-callable?]]
        [bcbio.variation.combine :only [combine-variants]]
        [bcbio.variation.config :only [load-config]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template get-vcf-source]])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]))

(defn- split-nocalls
  "Provide split VCFs of call and no-call variants for the given sample."
  [in-vcf sample ref out-dir]
  (letfn [(sample-only-vc [vc]
            (-> (VariantContextBuilder. vc)
                (.genotypes (GenotypesContext/create
                             (into-array [(-> vc .getGenotypes (.get sample))])))
                (.make)))
          (split-nocall-vc [vc]
            (when (empty? (:filters vc))
              (let [cur-vc (sample-only-vc (:vc vc))]
                [(if (.isNoCall (-> cur-vc .getGenotypes (.get sample))) :nocall :called)
                 cur-vc])))
          (set-header-to-sample [sample header]
            (VCFHeader. (.getMetaData header) (ordered-set sample)))]
    (let [out {:called (itx/add-file-part in-vcf (str sample "-called") out-dir)
               :nocall (itx/add-file-part in-vcf (str sample "-nocall") out-dir)}]
      (when (itx/needs-run? (vals out))
        (with-open [in-vcf-s (get-vcf-source in-vcf ref)]
          (write-vcf-w-template in-vcf out
                                (remove nil? (map split-nocall-vc (parse-vcf in-vcf-s)))
                                ref
                                :header-update-fn (partial set-header-to-sample sample))))
      out)))

(defn call-at-known-alleles
  "Do UnifiedGenotyper calling at known variation alleles."
  [site-vcf align-bam ref]
  (let [file-info {:out-vcf (itx/add-file-part site-vcf "wrefs")}
        args (concat ["-R" ref
                      "-o" :out-vcf
                      "-I" align-bam
                      "--alleles" site-vcf
                      "-L" site-vcf
                      "--genotyping_mode" "GENOTYPE_GIVEN_ALLELES"
                      "--output_mode" "EMIT_ALL_SITES"
                      "-stand_call_conf" "0.0"
                      "--genotype_likelihoods_model" "BOTH"]
                     (reduce #(concat %1 ["-A" %2]) [] std-annotations))]
    (broad/run-gatk "UnifiedGenotyper" args file-info {:out [:out-vcf]})
    (:out-vcf file-info)))

(defn recall-nocalls
  "Recall variations at no-calls in a sample using UnifiedGenotyper."
  [in-vcf sample align-bam ref & {:keys [out-dir]}]
  (let [out-file (itx/add-file-part in-vcf (str sample "-wrefs") out-dir)]
    (when (itx/needs-run? out-file)
      (let [{:keys [called nocall]} (split-nocalls in-vcf sample ref out-dir)
            ready-nocall (call-at-known-alleles nocall align-bam ref)]
        (fs/rename
         (combine-variants [called ready-nocall] ref :merge-type :full)
         out-file)))
    out-file))

(defn create-merged
  "Create merged VCF files with no-call/ref-calls for each of the inputs.
  Works at a higher level than `recall-nocalls` and does the work of
  preparing a set of all merged variants, then re-calling at non-missing positions."
  [vcfs align-bams vcf-configs ref & {:keys [out-dir intervals]}]
  (letfn [(merge-vcf [vcf sample all-vcf align-bam ref]
            (let [ready-vcf (combine-variants [vcf all-vcf] ref
                                              :merge-type :full :intervals intervals
                                              :out-dir out-dir :check-ploidy? false)]
              (recall-nocalls ready-vcf sample align-bam ref :out-dir out-dir)))]
    (let [merged (combine-variants vcfs ref :merge-type :minimal :intervals intervals
                                   :out-dir out-dir :check-ploidy? false)]
      (map (fn [[v b vcf-config]]
             (if (and (get vcf-config :refcalls false)
                      (not (nil? b)))
               (merge-vcf v (:name vcf-config) merged b ref)
               v))
           (map vector vcfs align-bams vcf-configs)))))

(defn convert-no-calls-w-callability
  "Convert no-calls into callable reference and real no-calls.
  Older functionality to re-call as reference when region is callable.
  Prefer `recall-nocalls`"
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

(defn -main [config-file]
  (let [config (load-config config-file)]
    (doseq [exp (config :experiments)]
      (doseq [call (exp :calls)]
        (recall-nocalls (:file call) (:name call) (:align call)
                        (:ref exp) :out-dir (get-in config [:dir :out]))))))