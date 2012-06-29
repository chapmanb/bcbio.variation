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
  (:use [clojure.java.io]
        [ordered.map :only [ordered-map]]
        [ordered.set :only [ordered-set]]
        [bcbio.variation.callable :only [get-callable-checker is-callable?]]
        [bcbio.variation.combine :only [combine-variants multiple-samples?]]
        [bcbio.variation.config :only [load-config]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template get-vcf-source
                                               get-vcf-header]])
  (:require [clojure.string :as string]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]))

(defn- split-nocalls
  "Provide split VCFs of call and no-call variants for the given sample."
  [in-vcf sample ref out-dir]
  (letfn [(sample-only-vc [vc]
            (-> (VariantContextBuilder. vc)
                (.genotypes (GenotypesContext/create
                             (into-array [(-> vc .getGenotypes (.get sample))])))
                (.attributes {})
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
        annotations ["DepthPerAlleleBySample"]
        args (concat ["-R" ref
                      "-o" :out-vcf
                      "-I" align-bam
                      "--alleles" site-vcf
                      "-L" site-vcf
                      "--genotyping_mode" "GENOTYPE_GIVEN_ALLELES"
                      "--output_mode" "EMIT_ALL_SITES"
                      "-stand_call_conf" "0.0"
                      "--genotype_likelihoods_model" "BOTH"]
                     (reduce #(concat %1 ["-A" %2]) [] annotations))]
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

(defn- split-vcf-sample-line
  "Split VCF line into shared attributes and sample specific genotypes.
  By default removes shared attributes which are no longer valid for split file."
  ([line remove-info-attrs?]
     (let [parts (string/split line #"\t")
           orig-shared (vec (take 9 parts))
           shared (if remove-info-attrs? (assoc orig-shared 7 ".") orig-shared)]
       (for [s (drop 9 parts)] (conj shared s))))
  ([line]
     (split-vcf-sample-line line true)))

(defn- split-vcf-to-samples-header
  "Split a multi-sample file to individual samples: writing the header."
  [vcf-iter out-files]
  (letfn [(not-chrom? [l] (not (.startsWith l "#CHROM")))]
    (let [std-header (string/join "\n" (take-while not-chrom? vcf-iter))]
      (doseq [[i xs] (map-indexed vector (-> (drop-while not-chrom? vcf-iter)
                                             first
                                             (split-vcf-sample-line false)))]
        (spit (get out-files i)
              (str std-header "\n" (string/join "\t" xs) "\n"))))))

(defn split-vcf-to-samples-variants
  "Split multi-sample file to individual samples: variant lines
  Avoids opening all output handles, instead writing to individual files.
  Blocks writes into groups to reduce opening file penalties."
  [vcf-iter out-files]
  (let [block-size 1000]
    (doseq [lines (partition-all block-size (drop-while #(.startsWith % "#") vcf-iter))]
      (let [sample-lines (reduce (fn [coll l]
                                   (reduce (fn [inner-coll [i xs]]
                                             (assoc inner-coll i (conj (get inner-coll i [])
                                                                       (string/join "\t" xs))))
                                           coll (map-indexed vector (split-vcf-sample-line l))))
                                 {} lines)]
        (doseq [[i xs] sample-lines]
          (spit (get out-files i)
                (str (string/join "\n" xs) "\n")
                :append true))))))

(defn split-vcf-to-samples
  "Create individual sample variant files from input VCF."
  [vcf-file & {:keys [out-dir]}]
  (let [samples (-> vcf-file get-vcf-header .getGenotypeSamples)
        out-files (into (ordered-map) (map (fn [x] [x (itx/add-file-part vcf-file x out-dir)])
                                           samples))]
    (when (itx/needs-run? (vals out-files))
      (with-open [rdr (reader vcf-file)]
        (itx/with-tx-files [tx-out-files out-files (keys out-files) []]
          (let [line-iter (line-seq rdr)]
            (split-vcf-to-samples-header line-iter (vec (vals tx-out-files)))
            (split-vcf-to-samples-variants line-iter (vec (vals tx-out-files)))))))
    out-files))

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
                       coll (split-vcf-to-samples vcf :out-dir out-dir)))
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
