(ns bcbio.variation.recall
  "Recall batched sets of variants containing no-call regions.
  Combined variant calls from batches contain regions called in
  some samples but not others. The approach:
    - Split sample into called and no-call variants
    - Re-call the no-call variants using the UnifiedGenotyper
    - Merge previously called and re-called into final set.
  http://www.broadinstitute.org/gsa/wiki/index.php/Merging_batched_call_sets"
  (:import [org.broadinstitute.sting.utils.variantcontext
            GenotypeBuilder VariantContextBuilder GenotypesContext]
           [org.broadinstitute.sting.utils.codecs.vcf VCFHeader])
  (:use [clojure.java.io]
        [ordered.map :only [ordered-map]]
        [ordered.set :only [ordered-set]]
        [bcbio.variation.callable :only [get-callable-checker is-callable?]]
        [bcbio.variation.combine :only [combine-variants fix-minimal-combined]]
        [bcbio.variation.config :only [load-config]]
        [bcbio.variation.filter.intervals :only [select-by-sample]]
        [bcbio.variation.haploid :only [diploid-calls-to-haploid]]
        [bcbio.variation.multisample :only [multiple-samples?]]
        [bcbio.variation.normalize :only [fix-vcf-sample remove-ref-alts]]
        [bcbio.variation.phasing :only [is-haploid?]])
  (:require [clojure.string :as string]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]
            [bcbio.variation.variantcontext :as gvc]))

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
            (VCFHeader. (.getMetaDataInInputOrder header) (ordered-set sample)))]
    (let [sample-str (if (.contains in-vcf sample) "" (str sample "-"))
          out {:called (itx/add-file-part in-vcf (str sample-str "called") out-dir)
               :nocall (itx/add-file-part in-vcf (str sample-str "nocall") out-dir)}]
      (when (itx/needs-run? (vals out))
        (with-open [in-vcf-iter (gvc/get-vcf-iterator in-vcf ref)]
          (gvc/write-vcf-w-template in-vcf out
                                    (remove nil? (map split-nocall-vc (gvc/parse-vcf in-vcf-iter)))
                                    ref
                                    :header-update-fn (partial set-header-to-sample sample))))
      out)))

(defn call-at-known-alleles
  "Do UnifiedGenotyper calling at known variation alleles."
  [site-vcf align-bam ref & {:keys [cores]}]
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
                      "-stand_emit_conf" "0.0"
                      "--max_deletion_fraction" "2.0"
                      "--min_indel_count_for_genotyping" "3"
                      "--genotype_likelihoods_model" "BOTH"]
                     (if cores ["-nt" (str cores)] [])
                     (reduce #(concat %1 ["-A" %2]) [] annotations))]
    (broad/index-bam align-bam)
    (broad/run-gatk "UnifiedGenotyper" args file-info {:out [:out-vcf]})
    (:out-vcf file-info)))

(defn recall-nocalls
  "Recall variations at no-calls in a sample using UnifiedGenotyper."
  [in-vcf sample call-name align-bam ref & {:keys [out-dir cores]}]
  (let [sample-str (if (.contains in-vcf call-name) "" (str call-name "-"))
        out-file (itx/add-file-part in-vcf (str sample-str "wrefs") out-dir)]
    (when (itx/needs-run? out-file)
      (let [{:keys [called nocall]} (split-nocalls in-vcf sample ref out-dir)
            prep-nocall (remove-ref-alts nocall ref)
            orig-nocall (call-at-known-alleles prep-nocall align-bam ref :cores cores)
            fix-nocall (fix-vcf-sample orig-nocall sample ref)
            ready-nocall (if (is-haploid? called ref)
                           (diploid-calls-to-haploid fix-nocall ref)
                           fix-nocall)
            combine-out (combine-variants [called ready-nocall] ref :merge-type :full
                                          :quiet-out? true)]
        (fs/rename combine-out out-file)
        (fs/rename (str combine-out ".idx") (str out-file ".idx"))))
    out-file))

(defn- no-recall-vcfs
  "Retrieve inputs VCFs not involved in preparing a recall VCF.
   This removes both the recalled input file, as well as
   individual inputs without recalling to avoid double counting."
  [all-vcfs vcf-configs]
  (let [remove-vcfs (set (map :file (filter :recall vcf-configs)))]
    (->> (interleave all-vcfs vcf-configs)
                         (partition 2)
                         (map (fn [[x c]]
                                (when-not (contains? remove-vcfs (:file c)) x)))
                         (remove nil?))))

;; ## Pick consensus variants

(defn- get-sample-call
  "Retrieve variant alleles for the sample, sorted in a stable order."
  [sample vc]
  (let [allele-order (->> (:alt-alleles vc)
                          (sort-by #(.getBaseString %))
                          (cons (:ref-allele vc))
                          (map-indexed vector)
                          (map reverse)
                          (map vec)
                          (into {}))]
    (->> (:genotypes vc)
         (filter #(= sample (:sample-name %)))
         first
         :alleles
         (sort-by allele-order))))

(defn- update-genotype-w-alleles
  "Update a genotype with the provided alleles."
  [vc sample alleles]
  (doto (-> vc :vc .getGenotypes GenotypesContext/copy)
    (.replace (-> (:vc vc)
                  .getGenotypes
                  (.get sample)
                  GenotypeBuilder.
                  (.alleles alleles)
                  (.phased false)
                  .make))))

(defn- update-vc-w-consensus
  "Update a variant context with consensus genotype from multiple inputs.
   Calculates the consensus set of calls, swapping calls to that if it
   exists. If there is no consensus default to the existing allele call."
  [vc sample input-vc-getter]
  (let [match-fn (juxt :start :ref-allele)
        max-alleles (->> (gvc/variants-in-region input-vc-getter vc)
                         (filter #(= (match-fn %) (match-fn vc)))
                         (cons vc)
                         (map (partial get-sample-call sample))
                         frequencies
                         (sort-by second >)
                         (partition-by second)
                         first)]
    (-> (VariantContextBuilder. (:vc vc))
        (.genotypes (update-genotype-w-alleles vc sample
                                               (if (= 1 (count max-alleles))
                                                 (ffirst max-alleles)
                                                 (get-sample-call sample vc))))
        .make)))

(defn- recall-w-consensus
  "Recall variants in a combined set of variants based on consensus of all inputs."
  [base-vcf input-vcfs sample ref-file]
  (let [out-file (itx/add-file-part base-vcf "consensus")]
    (when (itx/needs-run? out-file)
      (with-open [in-vcf-iter (gvc/get-vcf-iterator base-vcf ref-file)
                  input-vc-getter (apply gvc/get-vcf-retriever (cons ref-file input-vcfs))]
        (gvc/write-vcf-w-template base-vcf {:out out-file}
                                  (map #(update-vc-w-consensus % sample input-vc-getter)
                                       (gvc/parse-vcf in-vcf-iter))
                                  ref-file)))
    out-file))

(defn create-merged
  "Create merged VCF files with no-call/ref-calls for each of the inputs.
  Works at a higher level than `recall-nocalls` and does the work of
  preparing a set of all merged variants, then re-calling at non-missing positions."
  [vcfs align-bams exp & {:keys [out-dir intervals cores]}]
  (letfn [(merge-vcf [vcf call-name all-vcf align-bam ref]
            (-> [vcf all-vcf]
                (combine-variants ref :merge-type :full :intervals intervals
                                  :out-dir out-dir :check-ploidy? false)
                (recall-nocalls (:sample exp) call-name align-bam ref
                                :out-dir out-dir :cores cores)
                (recall-w-consensus (no-recall-vcfs vcfs (:calls exp))
                                    (:sample exp) ref)))]
    (let [min-merged (-> (combine-variants vcfs (:ref exp) :merge-type :minimal :intervals intervals
                                       :out-dir out-dir :check-ploidy? false
                                       :name-map (zipmap vcfs (map :name (:calls exp))))
                         (fix-minimal-combined vcfs (:ref exp)))]
      (map (fn [[v b vcf-config]]
             (if (and (get vcf-config :recall false)
                      (not (nil? b)))
               (let [merged (merge-vcf v (:name vcf-config) min-merged b (:ref exp))]
                 (if (get vcf-config :remove-refcalls true)
                   (select-by-sample (:sample exp) merged nil (:ref exp)
                                     :remove-refcalls true :ext "cleaned")
                   merged))
               v))
           (map vector vcfs align-bams (:calls exp))))))

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
  (let [samples (-> vcf-file gvc/get-vcf-header .getGenotypeSamples)
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
  (let [multi-files (filter multiple-samples? (set (map :file calls)))
        cur-samples (set (map :name calls))]
    (vals
     (reduce (fn [coll vcf]
               (reduce (fn [inner-coll [sample split-vcf]]
                         (if (contains? cur-samples sample)
                           (assoc-in inner-coll [sample :file] split-vcf)
                           inner-coll))
                       coll (split-vcf-to-samples vcf :out-dir out-dir)))
             (into (ordered-map) (map (fn [x] [(:name x) x]) calls))
             multi-files))))

(defn- remove-sample-info
  "Provide a cleaned VCF file without sample genotype information."
  [in-vcf out-dir]
  (letfn [(split-variant-line [line]
            (->> (string/split line #"\t")
                 (take 8)
                 (string/join "\t")))
          (process-line [line]
            (cond
             (.startsWith line "##fileformat") line
             (.startsWith line "##INFO") line
             (.startsWith line "#CHROM") (split-variant-line line)
             (.startsWith line "#") nil
             :else (split-variant-line line)))]
    (let [out-file (itx/add-file-part in-vcf "nosamples" out-dir)]
      (when (itx/needs-run? out-file)
        (with-open [rdr (reader in-vcf)
                    wtr (writer out-file)]
          (doseq [line (map process-line (line-seq rdr))]
            (when line
              (.write wtr (str line "\n"))))))
      out-file)))

(defn batch-combine-variants
  "Combine large numbers of variants via batches to avoid memory issues."
  [vcfs ref & {:keys [merge-type out-dir intervals unsafe name-map
                      base-ext check-ploidy? quiet-out? batch-size]
               :or {merge-type :unique
                    unsafe false
                    name-map {}
                    check-ploidy? true
                    batch-size 100}}]
  (letfn [(combine-w-args [xs]
            (combine-variants xs ref :merge-type merge-type :out-dir out-dir
                              :intervals intervals :unsafe unsafe :name-map name-map
                              :base-ext base-ext :check-ploidy? check-ploidy?
                              :quiet-out? quiet-out?))]
    (let [batch-vcfs (map combine-w-args (partition-all batch-size vcfs))]
      (combine-w-args batch-vcfs))))

(defn- do-recall-exp
  "Perform recalling on all specific inputs in an experiment"
  [exp out-dir config]
  (let [recall-vcfs (map (fn [call]
                           (recall-nocalls (:file call) (:sample exp) (:name call) (:align call)
                                           (:ref exp) :out-dir out-dir
                                           :cores (get-in config [:resources :cores])))
                         (split-config-multi (:calls exp) (:ref exp) out-dir))
        clean-multi (map #(remove-sample-info % out-dir)
                         (filter multiple-samples? (set (map :file (:calls exp)))))]
    (batch-combine-variants (concat clean-multi recall-vcfs) (:ref exp) :merge-type :full
                            :quiet-out? true :check-ploidy? false)))

(defn convert-no-calls-w-callability
  "Convert no-calls into callable reference and real no-calls.
  Older functionality to re-call as reference when region is callable.
  Prefer `recall-nocalls`"
  [in-vcf align-bam ref & {:keys [out-dir intervals num-alleles]}]
  (letfn [(ref-genotype [g vc]
            (doto (-> vc :vc .getGenotypes GenotypesContext/copy)
              (.replace
               (-> (GenotypeBuilder. (:genotype g))
                   (.alleles (repeat (if (nil? num-alleles)
                                       (count (:alleles g))
                                       num-alleles)
                                     (:ref-allele vc)))
                   .make))))
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
            (for [vc (gvc/parse-vcf vcf-source)]
              [:out (maybe-callable-vc vc call-source)]))]
    (let [out-file (itx/add-file-part in-vcf "wrefs")]
      (when (itx/needs-run? out-file)
        (with-open [in-vcf-iter (gvc/get-vcf-iterator in-vcf ref)
                    call-source (get-callable-checker align-bam ref :out-dir out-dir
                                                      :intervals intervals)]
          (gvc/write-vcf-w-template in-vcf {:out out-file}
                                    (convert-vcs in-vcf-iter call-source) ref)))
      out-file)))

(defn -main [config-file]
  (let [config (load-config config-file)
        out-dir (get-in config [:dir :out])]
    (doseq [exp (:experiments config)]
      (do-recall-exp exp out-dir config))))
