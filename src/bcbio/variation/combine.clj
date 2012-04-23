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
        [bcbio.variation.normalize :only [prep-vcf clean-problem-vcf]])
  (:require [fs.core :as fs]
            [clojure.string :as string]
            [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]))

(defn gatk-cl-intersect-intervals
  "Supply GATK commandline arguments for interval files, merging via intersection."
  [intervals]
  (cond
   (nil? intervals) []
   (coll? intervals) (concat (flatten (map #(list "-L" %) intervals))
                             ["--interval_set_rule" "INTERSECTION"])
   :else ["-L", intervals]))

(defn combine-variants
  "Combine multiple variant files with GATK CombineVariants."
  [vcfs ref & {:keys [merge-type out-dir intervals unsafe name-map base-ext]
               :or {merge-type :unique unsafe false name-map {}}}]
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
                       (flatten (map #(list (str "--variant:" (unique-name %)) %) vcfs))
                       (gatk-cl-intersect-intervals intervals)
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
  [in-vcf align-bam ref & {:keys [out-dir] :or {out-dir nil}}]
  (let [out-file (itx/add-file-part in-vcf "wrefs")
        [is-callable? call-source] (callable-checker align-bam ref :out-dir out-dir)]
    (letfn [(ref-genotype [g vc]
              (doto (-> vc :vc .getGenotypes GenotypesContext/copy)
                (.replace
                 (Genotype/modifyAlleles (:genotype g)
                                         (repeat (count (:alleles g))
                                                 (:ref-allele vc))))))
            (maybe-callable-vc [vc]
              {:pre (= 1 (count (:genotypes vc)))}
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
      (if (itx/needs-run? out-file)
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
    (let [vcf-samples (with-open [vcf-source (get-vcf-source in-vcf ref-file)]
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
                     (if-not (nil? intervals) ["-L" intervals] []))]
    (if-not (fs/exists? base-dir)
      (fs/mkdirs base-dir))
    (broad/run-gatk "SelectVariants" args file-info {:out [:out-vcf]})
    (:out-vcf file-info)))

(defn create-merged
  "Create merged VCF files with no-call/ref-calls for each of the inputs."
  [vcfs align-bams do-merges ref & {:keys [out-dir intervals]
                                    :or {out-dir nil intervals nil}}]
  (letfn [(merge-vcf [vcf all-vcf align-bam ref]
            (let [ready-vcf (combine-variants [vcf all-vcf] ref
                                              :merge-type :full :intervals intervals
                                              :out-dir out-dir)]
              (convert-no-calls ready-vcf align-bam ref :out-dir out-dir)))]
    (let [merged (combine-variants vcfs ref :merge-type :minimal :intervals intervals
                                   :out-dir out-dir)]
      (map (fn [[v b merge?]] (if merge? (merge-vcf v merged b ref) v))
           (map vector vcfs align-bams do-merges)))))

(defn- dirty-prep-work
  "Prepare input file for comparisons based on configuration:
    - Selecting a single sample from multi-sample files
    - Resorting and fixing chromosome naming
    - Removing reference call genotypes
   This organizes the logic which get convoluted for different cases.
   The approach is to select a single sample and remove refcalls if we have
   a multiple sample file, so the sample name will be correct."
  [in-file call exp out-dir out-fname]
  (letfn [(run-sample-select [in-file]
            (select-by-sample (:sample exp) in-file (:name call)
                              (get call :ref (:ref exp))
                              :out-dir out-dir
                              :remove-refcalls (get call :remove-refcalls false)))]
    (let [sample-file (if (multiple-samples? in-file)
                        (run-sample-select in-file)
                        in-file)
          prep-file (if (true? (:prep call))
                      (prep-vcf sample-file (:ref exp) (:sample exp) :out-dir out-dir
                                :out-fname out-fname :orig-ref-file (:ref call))
                      sample-file)
          noref-file (if (and (not (multiple-samples? in-file)) (:remove-refcalls call))
                       (run-sample-select prep-file)
                       prep-file)]
      noref-file)))

(defn gatk-normalize
  "Prepare call information for VCF comparisons by normalizing through GATK.
  Handles:

   1. Combining multiple input files
   2. Fixing reference and sample information.
   3. Splitting combined MNPs into phased SNPs"
  [call exp out-dir]
  (if-not (fs/exists? out-dir)
    (fs/mkdirs out-dir))
  (letfn [(merge-call-files [call in-files]
            (combine-variants in-files (get call :ref (:ref exp))
                              :merge-type :full :out-dir out-dir
                              :unsafe true))]
    (let [out-fname (format "%s-%s.vcf" (:sample exp) (:name call))
          in-files (if (coll? (:file call)) (:file call) [(:file call)])
          clean-files (map #(if-not (:preclean call) %
                                    (clean-problem-vcf % :out-dir out-dir))
                           in-files)
          merge-file (if (> (count clean-files) 1)
                       (merge-call-files call clean-files)
                       (first clean-files))
          prep-file (dirty-prep-work merge-file call exp out-dir out-fname)]
      (assoc call :file (if (true? (get call :normalize true))
                          (normalize-variants prep-file (:ref exp) out-dir
                                              :out-fname out-fname)
                          prep-file)))))
