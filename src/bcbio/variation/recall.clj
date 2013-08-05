(ns bcbio.variation.recall
  "Recall batched sets of variants using consensus from multiple inputs."
  (:import [org.broadinstitute.variant.variantcontext VariantContextBuilder]
           [org.broadinstitute.variant.vcf VCFHeader])
  (:use [clojure.java.io]
        [ordered.map :only [ordered-map]]
        [ordered.set :only [ordered-set]]
        [bcbio.variation.callable :only [get-callable-checker is-callable?]]
        [bcbio.variation.combine :only [combine-variants fix-minimal-combined]]
        [bcbio.variation.filter.intervals :only [select-by-sample]]
        [bcbio.variation.multisample :only [multiple-samples?]])
  (:require [clojure.string :as string]
            [bcbio.run.itx :as itx]
            [bcbio.variation.filter.attr :as attr]
            [bcbio.variation.variantcontext :as gvc]))

;; ## Utilities

(defn- set-header-to-samples [samples _ header]
  (VCFHeader. (.getMetaDataInInputOrder header) samples))

(defn- no-recall-vcfs
  "Retrieve inputs VCFs not involved in preparing a recall VCF.
   Avoid double pulling inputs with the same initial call files."
  [all-vcfs vcf-configs]
  (let [include-names (->> vcf-configs
                           (group-by :file)
                           (map second)
                           (map first)
                           (map :name)
                           set)]
    (->> (interleave all-vcfs vcf-configs)
         (partition 2)
         (map (fn [[x c]] (when (contains? include-names (:name c)) x)))
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
                          (into {}))
        g (->> (:genotypes vc)
                   (filter #(= sample (:sample-name %)))
                   first)]
    {:sample-name sample
     :qual (:qual vc)
     :vc-type (:type vc)
     :call-type (:type g)
     :ref-allele (:ref-allele vc)
     :alleles (sort-by allele-order (:alleles g))
     :attributes (select-keys (:attributes g) ["PL" "DP" "AD" "PVAL"])
     :has-likelihood (if (seq (get-in g [:attributes "PL"])) 1 0)
     :attr-count (+ (if (seq (get-in g [:attributes "PL"])) 1 0)
                    (if (seq (get-in g [:attributes "PVAL"])) 1 0)
                    (if (seq (get-in g [:attributes "AD"])) 1 0)
                    (if (get-in g [:attributes "DP"]) 1 0))
     :pl nil ;(attr/get-vc-attr vc "PL" nil)
     }))

(defn- best-supported-alleles
  "Retrieve alleles with best support from multiple inputs.
   Use posterior likelihoods and quality scores to rank results
   with the same alleles and counts. We rank by total number of
   calls identified. We break ties in favor of homozygous calls
   if there isn't a consensus on het/hom variant calls, since
   we've failed to establish power for calling more difficult hets."
  [alleles]
  (letfn [(safe-sum [xs k]
            (apply + (remove nil? (map k xs))))
          (sum-plus-call-type [i xs]
            (let [pls (safe-sum xs :pl)
                  represent-x (last (sort-by #(vector (:has-likelihood %)
                                                      (or (:pl %) (- Integer/MIN_VALUE))
                                                      (:attr-count %))
                                             xs))
                  call-code (if (= "HET" (:call-type represent-x)) 0 1)]
              [(count xs) call-code (- pls) i represent-x]))]
    (->> alleles
         (group-by :alleles)
         (map second)
         (map-indexed sum-plus-call-type)
         sort
         last ; Best item
         last ; Extract the alleles
         )))

(defn- best-supported-sample-gs
  "Identify the best supported genotypes for a sample from multiple variant calls."
  [vcs sample]
  (->> vcs
       (map (partial get-sample-call sample))
       best-supported-alleles))

(defn- update-vc-w-consensus
  "Update a variant context with consensus genotype from multiple inputs.
   Calculates the consensus set of calls, swapping calls to that if it
   exists. If there is no consensus default to the existing allele call."
  [vc samples input-vc-getter]
  (let [match-fn (juxt :start :ref-allele)
        other-vcs (filter #(= (match-fn %) (match-fn vc))
                          (gvc/variants-in-region input-vc-getter vc))
        most-likely-gs (->> samples
                            (map (partial best-supported-sample-gs other-vcs))
                            (remove nil?))]
    (when (seq most-likely-gs)
      (-> (VariantContextBuilder. (:vc vc))
          (.alleles (conj (set (remove #(.isNoCall %) (mapcat :alleles most-likely-gs)))
                          (:ref-allele vc)))
          (.genotypes (gvc/create-genotypes most-likely-gs :attrs #{"PL" "PVAL" "DP" "AD"}))
          .make))))

(defn- recall-w-consensus
  "Recall variants in a combined set of variants based on consensus of all inputs."
  [base-vcf input-vcfs sample ref-file]
  (let [out-file (itx/add-file-part base-vcf "consensus")]
    (when (itx/needs-run? out-file)
      (with-open [in-vcf-iter (gvc/get-vcf-iterator base-vcf ref-file)
                  input-vc-getter (apply gvc/get-vcf-retriever (cons ref-file input-vcfs))]
        (let [samples (into (ordered-set) (if sample
                                            [sample]
                                            (-> input-vcfs first gvc/get-vcf-header .getGenotypeSamples)))]
          (gvc/write-vcf-w-template base-vcf {:out out-file}
                                    (remove nil? (map #(update-vc-w-consensus % samples input-vc-getter)
                                                      (gvc/parse-vcf in-vcf-iter)))
                                    ref-file
                                    :header-update-fn (partial set-header-to-samples samples)))))
    out-file))

(defn- get-min-merged
  "Retrieve a minimal merged file with calls from input VCFs."
  [vcfs exp out-dir intervals]
  (-> (combine-variants vcfs (:ref exp) :merge-type :minimal :intervals intervals
                        :out-dir out-dir :check-ploidy? false
                        :name-map (zipmap vcfs (map :name (:calls exp))))
      (fix-minimal-combined vcfs (:ref exp))))

(defmulti recall-vcf
  "Recall missing calls, handling merging or consensus based approaches"
  (fn [in-info & _]
    (if (some nil? [in-info (:bam in-info) (:file in-info)])
      :consensus
      (keyword (get in-info :approach :consensus)))))

(defmethod recall-vcf :consensus
  ^{:doc "Provide recalling of nocalls based on consensus from all inputs."}
  [in-info vcfs exp out-dir intervals]
  (-> vcfs
      (get-min-merged exp out-dir intervals)
      (recall-w-consensus (no-recall-vcfs vcfs (:calls exp))
                          (:sample exp) (:ref exp))))

(defn create-merged
  "Create merged VCF files with no-call/ref-calls for each of the inputs.
  Works at a higher level than `recall-vcf` and does the work of
  preparing a set of all merged variants, then re-calling at non-missing positions."
  [vcfs align-bams exp & {:keys [out-dir intervals cores]}]
  (map (fn [[v b vcf-config]]
         (if (get vcf-config :recall false)
           (let [base-info {:name (:name vcf-config)
                            :approach (get-in exp [:params :recall-approach] :consensus)
                            :file v :bam b}
                 merged (recall-vcf base-info vcfs exp out-dir intervals)]
             (if (and (:sample exp) (get vcf-config :remove-refcalls true))
               (select-by-sample (:sample exp) merged nil (:ref exp)
                                 :remove-refcalls true :ext "cleaned")
               merged))
           v))
       (map vector vcfs align-bams (:calls exp))))

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

(defn convert-no-calls-w-callability
  "Convert no-calls into callable reference and real no-calls.
  Older functionality to re-call as reference when region is callable.
  Prefer `recall-vcf`"
  [in-vcf align-bam ref & {:keys [out-dir intervals num-alleles]}]
  (letfn [(maybe-callable-vc [vc call-source]
            {:pre (= 1 (:num-samples vc))}
            (let [g (-> vc :genotypes first)]
              (if (.isNoCall (-> g :alleles first))
                (if (is-callable? call-source (:chr vc) (:start vc) (:end vc))
                  (gvc/genotypes->refcall vc)
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
