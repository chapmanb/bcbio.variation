
(ns bcbio.variation.combine
  "Combine variant files, handling no-calls versus reference calls

   1. Combine the variants to create a merged set of positions to call at
   2. For each variant file:
      a. Generate callability at each position
      b. Combine original calls with merged positions
      c. Walk through each no-call and set as reference if callable"
  (:import [org.broadinstitute.sting.utils.variantcontext
            Genotype VariantContextBuilder GenotypesContext])
  (:use [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template]]
        [bcbio.variation.callable :only [callable-checker]]
        [bcbio.variation.complex :only [normalize-variants]]
        [bcbio.variation.normalize :only [prep-vcf]])
  (:require [fs.core :as fs]
            [clojure.string :as string]
            [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]))

(defn combine-variants
  "Combine two variant files with GATK CombineVariants."
  [vcfs ref & {:keys [merge-type out-dir intervals]
               :or {merge-type :unique out-dir nil intervals nil}}]
  (letfn [(unique-name [f]
            (-> f fs/base-name itx/file-root))]
    (let [base-dir (if (nil? out-dir) (fs/parent (first vcfs)) out-dir)
          file-info {:out-vcf (str (fs/file base-dir
                                            (itx/add-file-part (-> vcfs
                                                                   first
                                                                   fs/base-name
                                                                   itx/remove-zip-ext)
                                                               (case merge-type
                                                                     :minimal "mincombine"
                                                                     :full "fullcombine"
                                                                     "combine"))))}
          args (concat ["-R" ref
                        "-o" :out-vcf
                        "--rod_priority_list" (string/join "," (map unique-name vcfs))]
                       (flatten (map #(list (str "--variant:" (unique-name %)) %) vcfs))
                       (cond
                        (nil? intervals) []
                        (coll? intervals) (concat (flatten (map #(list "-L" %) intervals))
                                                  ["--interval_set_rule" "INTERSECTION"])
                        :else ["-L", intervals])
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
        is-callable? (callable-checker align-bam ref :out-dir out-dir)]
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
            (convert-vcs [in-file]
              (for [vc (parse-vcf in-file)]
                [:out (maybe-callable-vc vc)]))]
      (if (itx/needs-run? out-file)
        (write-vcf-w-template in-vcf {:out out-file} (convert-vcs in-vcf) ref))
      out-file)))

(defn select-by-sample
  "Select only the sample of interest from input VCF files."
  [sample call ref & {:keys [out-dir intervals]
                      :or {out-dir nil intervals nil}}]
  (let [base-dir (if (nil? out-dir) (fs/parent (:file call)) out-dir)
        file-info {:out-vcf (str (fs/file base-dir
                                          (format "%s-%s.vcf" sample (:name call))))}
        args (concat ["-R" ref
                      "--sample_name" sample
                      "--variant" (:file call)
                      "--out" :out-vcf]
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
                                              :merge-type :full :intervals intervals)]
              (convert-no-calls ready-vcf align-bam ref :out-dir out-dir)))]
    (let [merged (combine-variants vcfs ref :merge-type :minimal :intervals intervals)]
      (map (fn [[v b merge?]] (if merge? (merge-vcf v merged b ref) v))
           (map vector vcfs align-bams do-merges)))))

(defn gatk-normalize
  "Prepare call information for VCF comparisons by normalizing through GATK.
  Handles:

   1. Combining multiple input files
   2. Fixing reference and sample information.
   3. Splitting combined MNPs into phased SNPs"
  [call exp out-dir]
  (if-not (fs/exists? out-dir)
    (fs/mkdirs out-dir))
  (letfn [(merge-call-files [call]
            (combine-variants (:file call) (:ref exp)
                              :merge-type :full :out-dir out-dir))]
    (let [merge-file (if (coll? (:file call))
                       (merge-call-files call)
                       (:file call))
          prep-file (if (true? (:prep call))
                      (prep-vcf merge-file (:ref exp) (:sample exp) :out-dir out-dir)
                      merge-file)]
      (assoc call :file (normalize-variants prep-file (:ref exp) out-dir)))))
