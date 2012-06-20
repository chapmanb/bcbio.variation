(ns bcbio.variation.combine
  "Combine variant files, handling no-calls versus reference calls

   1. Combine the variants to create a merged set of positions to call at
   2. For each variant file:
      a. Generate callability at each position
      b. Combine original calls with merged positions
      c. Walk through each no-call and set as reference if callable"
  (:import [org.broadinstitute.sting.utils.variantcontext
            Genotype VariantContextBuilder GenotypesContext])
  (:use [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template get-vcf-source
                                               get-vcf-header]]
        [bcbio.variation.callable :only [callable-checker]]
        [bcbio.variation.complex :only [normalize-variants]]
        [bcbio.variation.haploid :only [diploid-calls-to-haploid]]
        [bcbio.variation.normalize :only [prep-vcf clean-problem-vcf]]
        [bcbio.variation.phasing :only [is-haploid?]])
  (:require [fs.core :as fs]
            [clojure.string :as string]
            [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]))

(defn combine-variants
  "Combine multiple variant files with GATK CombineVariants.
   Only correctly handles all-by-all comparisons with the same ploidy level."
  [vcfs ref & {:keys [merge-type out-dir intervals unsafe name-map base-ext check-ploidy? quiet-out?]
               :or {merge-type :unique
                    unsafe false
                    name-map {}
                    check-ploidy? true}}]
  (when (and check-ploidy?
             (> (count (set (remove nil? (map #(is-haploid? % ref) vcfs)))) 1))
    (throw (Exception. (format "Haploid and non-haploid combinations not supported: %s %s"
                               (vec vcfs) (vec (map #(is-haploid? % ref) vcfs))))))
  (letfn [(unique-name [f]
            (get name-map f
                 (-> f fs/base-name itx/file-root)))]
    (let [base-dir (if (nil? out-dir) (fs/parent (first vcfs)) out-dir)
          full-base-name (-> vcfs first fs/base-name itx/remove-zip-ext)
          base-name (if (nil? base-ext) full-base-name
                        (format "%s-%s.vcf" (first (string/split full-base-name #"-"))
                                base-ext))
          file-info {:out-vcf (str (fs/file base-dir
                                            (itx/add-file-part base-name
                                                               (case merge-type
                                                                     :minimal "mincombine"
                                                                     :full "fullcombine"
                                                                     "combine"))))}
          args (concat ["-R" ref
                        "-o" :out-vcf
                        "--rod_priority_list" (string/join "," (map unique-name vcfs))]
                       (if unsafe ["--unsafe" "ALLOW_SEQ_DICT_INCOMPATIBILITY"] [])
                       (if quiet-out? ["--suppressCommandLineHeader" "--setKey" "null"] [])
                       (flatten (map #(list (str "--variant:" (unique-name %)) %) vcfs))
                       (broad/gatk-cl-intersect-intervals intervals ref)
                       (case merge-type
                             :full ["--genotypemergeoption" "PRIORITIZE"]
                             :unique ["--genotypemergeoption" "UNIQUIFY"]
                             :minimal ["--sites_only" "--minimalVCF"]))]
      (if-not (fs/exists? base-dir)
        (fs/mkdirs base-dir))
      (broad/run-gatk "CombineVariants" args file-info {:out [:out-vcf]})
      (:out-vcf file-info))))

(defn convert-no-calls
  "Convert no-calls into callable reference and real no-calls."
  [in-vcf align-bam ref & {:keys [out-dir intervals num-alleles]}]
  (let [out-file (itx/add-file-part in-vcf "wrefs")
        [is-callable? call-source] (callable-checker align-bam ref :out-dir out-dir
                                                     :intervals intervals)]
    (letfn [(ref-genotype [g vc]
              (doto (-> vc :vc .getGenotypes GenotypesContext/copy)
                (.replace
                 (Genotype/modifyAlleles (:genotype g)
                                         (repeat (if (nil? num-alleles)
                                                   (count (:alleles g))
                                                   num-alleles)
                                                 (:ref-allele vc))))))
            (maybe-callable-vc [vc]
              {:pre (= 1 (:num-samples vc))}
              (let [g (-> vc :genotypes first)]
                (if (.isNoCall (-> g :alleles first))
                  (if (is-callable? (:chr vc) (:start vc) (:end vc))
                    (-> (VariantContextBuilder. (:vc vc))
                        (.genotypes (ref-genotype g vc))
                        (.make))
                    (-> (VariantContextBuilder. (:vc vc))
                        (.filters #{"NotCallable"})
                        (.make)))
                  (:vc vc))))
            (convert-vcs [vcf-source]
              (for [vc (parse-vcf vcf-source)]
                [:out (maybe-callable-vc vc)]))]
      (when (itx/needs-run? out-file)
        (with-open [in-vcf-s (get-vcf-source in-vcf ref)
                    _ call-source]
          (write-vcf-w-template in-vcf {:out out-file} (convert-vcs in-vcf-s) ref)))
      out-file)))

(defn multiple-samples?
  "Check if the input VCF file has multiple genotyped samples."
  [in-file & {:keys [sample]}]
  (let [samples (-> in-file get-vcf-header .getGenotypeSamples)]
    (or (> (count samples) 1)
        (and (not (nil? sample))
             (not (contains? (set samples) sample))))))

(defn vcf-sample-name
  "Retrieve the sample name in a provided VCF file, allowing for partial matches."
  [sample in-vcf ref-file]
  (letfn [(sample-match [x choices]
            (let [do-match (filter #(when (.contains % x) %) choices)]
              (when (= 1 (count do-match))
                (first do-match))))]
    (let [vcf-samples (with-open [vcf-source (get-vcf-source in-vcf ref-file
                                                             :ensure-safe false)]
                        (-> vcf-source .getHeader .getGenotypeSamples set))]
      (if (contains? vcf-samples sample)
        sample
        (sample-match sample vcf-samples)))))

(defn select-by-sample
  "Select only the sample of interest from input VCF files."
  [sample in-file name ref & {:keys [out-dir intervals remove-refcalls]
                              :or {remove-refcalls false}}]
  (let [base-dir (if (nil? out-dir) (fs/parent in-file) out-dir)
        file-info {:out-vcf (str (fs/file base-dir
                                          (format "%s-%s.vcf" sample name)))}
        args (concat ["-R" ref
                      "--sample_name" (vcf-sample-name sample in-file ref)
                      "--variant" in-file
                      "--unsafe" "ALLOW_SEQ_DICT_INCOMPATIBILITY"
                      "--out" :out-vcf]
                     (if remove-refcalls ["--excludeNonVariants" "--excludeFiltered"] [])
                     (broad/gatk-cl-intersect-intervals intervals ref))]
    (if-not (fs/exists? base-dir)
      (fs/mkdirs base-dir))
    (broad/run-gatk "SelectVariants" args file-info {:out [:out-vcf]})
    (:out-vcf file-info)))

(defn create-merged
  "Create merged VCF files with no-call/ref-calls for each of the inputs."
  [vcfs align-bams do-merges ref & {:keys [out-dir intervals]}]
  (letfn [(merge-vcf [vcf all-vcf align-bam ref]
            (let [ready-vcf (combine-variants [vcf all-vcf] ref
                                              :merge-type :full :intervals intervals
                                              :out-dir out-dir :check-ploidy? false)
                  num-alleles (when (is-haploid? vcf ref) 1)]
              (convert-no-calls ready-vcf align-bam ref :out-dir out-dir
                                :intervals intervals :num-alleles num-alleles)))]
    (let [merged (combine-variants vcfs ref :merge-type :minimal :intervals intervals
                                   :out-dir out-dir :check-ploidy? false)]
      (map (fn [[v b merge?]] (if merge? (merge-vcf v merged b ref) v))
           (map vector vcfs align-bams do-merges)))))

(defn- genome-safe-intervals
  "Check if interval BED files overlap with current analysis genome build.
  This is useful when an input VCF is from an alternate genome and needs
  conversion. In this case we shouldn't yet be using interval selection."
  [intervals ref exp]
  (if (or (nil? ref) (= ref (:ref exp)))
    intervals
    []))

(defn- dirty-prep-work
  "Prepare input file for comparisons based on configuration:
    - Selecting a single sample from multi-sample files
    - Resorting and fixing chromosome naming
    - Removing reference call genotypes
   This organizes the logic which get convoluted for different cases.
   The approach is to select a single sample and remove refcalls if we have
   a multiple sample file, so the sample name will be correct."
  [in-file call exp intervals out-dir out-fname]
  (letfn [(run-sample-select [in-file ref ext]
            (select-by-sample (:sample exp) in-file (str (:name call) ext)
                              ref :out-dir out-dir
                              :intervals (genome-safe-intervals intervals ref exp)
                              :remove-refcalls (get call :remove-refcalls false)))]
    (let [sample-file (if (multiple-samples? in-file)
                        (run-sample-select in-file (get call :ref (:ref exp)) "")
                        in-file)
          prep-file (if (true? (:prep call))
                      (prep-vcf sample-file (:ref exp) (:sample exp) :out-dir out-dir
                                :out-fname out-fname :orig-ref-file (:ref call)
                                :config call)
                      sample-file)
          hap-file (if (true? (:make-haploid call))
                     (diploid-calls-to-haploid prep-file (:ref exp) :out-dir out-dir)
                     prep-file)
          noref-file (if (or (and (not (multiple-samples? in-file)) (:remove-refcalls call))
                             (and (not (nil? (:ref call))) (not (empty? intervals))))
                       (run-sample-select hap-file (:ref exp) "-noref")
                       hap-file)]
      noref-file)))

(defn gatk-normalize
  "Prepare call information for VCF comparisons by normalizing through GATK.
  Handles:

   1. Combining multiple input files
   2. Fixing reference and sample information.
   3. Splitting combined MNPs into phased SNPs"
  [call exp intervals out-dir transition]
  (if-not (fs/exists? out-dir)
    (fs/mkdirs out-dir))
  (letfn [(merge-call-files [call in-files]
            (let [ref (get call :ref (:ref exp))]
              (combine-variants in-files ref
                                :merge-type :full :out-dir out-dir
                                :intervals (genome-safe-intervals intervals ref exp)
                                :unsafe true)))]
    (let [out-fname (format "%s-%s.vcf" (:sample exp) (:name call))
          in-files (if (coll? (:file call)) (:file call) [(:file call)])
          _ (transition :clean (str "Cleaning input VCF: " (:name call)))
          clean-files (vec (map #(if-not (:preclean call) %
                                         (clean-problem-vcf % :out-dir out-dir))
                                in-files))
          _ (transition :merge (str "Merging multiple input files: " (:name call)))
          merge-file (if (> (count clean-files) 1)
                       (merge-call-files call clean-files)
                       (first clean-files))
          _ (transition :prep
                        (str "Resorting to comparison reference and selecting samples: "
                             (:name call)))
          prep-file (dirty-prep-work merge-file call exp intervals out-dir out-fname)]
      (transition :normalize (str "Normalize MNP and indel variants: " (:name call)))
      (assoc call :file (if (true? (get call :normalize true))
                          (normalize-variants prep-file (:ref exp) out-dir
                                              :out-fname out-fname)
                          prep-file)))))
