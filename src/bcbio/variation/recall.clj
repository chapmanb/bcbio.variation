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
  (:use [ordered.map :only [ordered-map]]
        [ordered.set :only [ordered-set]]
        [bcbio.variation.annotation :only [std-annotations]]
        [bcbio.variation.callable :only [get-callable-checker is-callable?]]
        [bcbio.variation.combine :only [combine-variants multiple-samples?]]
        [bcbio.variation.config :only [load-config]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template get-vcf-source
                                               get-vcf-header]])
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
          (set-header-to-sample [sample _ header]
            (VCFHeader. (.getMetaData header) (ordered-set sample)))]
    (let [sample-str (if (.contains in-vcf sample) "" (str sample "-"))
          out {:called (itx/add-file-part in-vcf (str sample-str "called") out-dir)
               :nocall (itx/add-file-part in-vcf (str sample-str "nocall") out-dir)}]
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
  (let [sample-str (if (.contains in-vcf sample) "" (str sample "-"))
        out-file (itx/add-file-part in-vcf (str sample-str "wrefs") out-dir)]
    (when (itx/needs-run? out-file)
      (let [{:keys [called nocall]} (split-nocalls in-vcf sample ref out-dir)
            ready-nocall (call-at-known-alleles nocall align-bam ref)]
        (fs/rename
         (combine-variants [called ready-nocall] ref :merge-type :full :quiet-out? true)
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

(defn- split-vcf-to-samples-batch
  "Create batch of individual sample variant files from input VCF"
  [vcf-file samples ref-file out-dir]
  (letfn [(get-one-sample-vcs [vc]
            (map (fn [g] [(:sample-name g)
                          (-> (VariantContextBuilder. (:vc vc))
                              (.genotypes (GenotypesContext/create (into-array [(:genotype g)])))
                              .make)])
                 (filter #(contains? samples (:sample-name %)) (:genotypes vc))))
          (set-header-to-sample [sample header]
            (VCFHeader. (.getMetaData header) (ordered-set sample)))]
    (let [out-files (reduce (fn [coll x]
                              (assoc coll x (itx/add-file-part vcf-file x out-dir)))
                            {} samples)]
      (when (itx/needs-run? (vals out-files))
        (with-open [vcf-source (get-vcf-source vcf-file ref-file)]
          (write-vcf-w-template vcf-file out-files
                                (partition 2 (flatten (map get-one-sample-vcs (parse-vcf vcf-source))))
                                ref-file
                                :header-update-fn set-header-to-sample)))
      out-files)))

(defn split-vcf-to-samples
  "Create individual sample variant files from input VCF.
  Handles batching inputs into groups to avoid too-many-file-open errors"
  [vcf-file ref-file & {:keys [out-dir max-files]
                        :or {max-files 100}}]
  (reduce (fn [out samples]
            (merge out (split-vcf-to-samples-batch vcf-file (set samples) ref-file out-dir)))
          {}
          (partition-all max-files (-> vcf-file get-vcf-header .getGenotypeSamples))))

(defn- split-config-multi
  "Split multiple sample inputs into individual samples before processing.
  This helps reduce the load on selecting from huge multi-sample files.
  Returns a list of configured calls with multi-samples set to individually
  separated input files."
  [calls ref out-dir]
  (let [multi-files (filter multiple-samples? (set (map :file calls)))]
    (vals
     (reduce (fn [coll vcf]
               (reduce (fn [inner-coll [sample split-vcf]]
                         (assoc-in inner-coll [sample :file] split-vcf))
                       coll (split-vcf-to-samples vcf ref :out-dir out-dir)))
             (into (ordered-map) (map (fn [x] [(:name x) x]) calls))
             multi-files))))

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
  (let [config (load-config config-file)
        out-dir (get-in config [:dir :out])
        recall-vcfs (pmap (fn [[exp call]]
                            (recall-nocalls (:file call) (:name call) (:align call)
                                            (:ref exp) :out-dir out-dir))
                          (apply concat
                                 (for [exp (:experiments config)]
                                   (for [call (split-config-multi (:calls exp) (:ref exp) out-dir)]
                                     [exp call]))))]
    (println recall-vcfs)))
