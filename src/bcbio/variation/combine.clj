;; Combine variant files, handling no-calls versus reference calls
;; 1. Combine the variants to create a merged set of positions to call at
;; 2. For each variant file:
;;    a. Generate callability at each position
;;    b. Combine original calls with merged positions
;;    c. Walk through each no-call and set as reference if callable

(ns bcbio.variation.combine
  (:import [org.broadinstitute.sting.utils.variantcontext
            Genotype VariantContextBuilder GenotypesContext])
  (:use [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template]]
        [bcbio.variation.callable :only [callable-checker]]
        [bcbio.variation.complex :only [normalize-variants]]
        [clojure.string :only [join]])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]))

(defn combine-variants [vcfs ref & {:keys [merge-type out-dir intervals]
                                    :or {merge-type :unique
                                         out-dir nil
                                         intervals nil}}]
  "Combine two variant files with GATK CombineVariants."
  (letfn [(unique-name [f]
            (-> f fs/base-name itx/file-root))]
    (let [base-dir (if (nil? out-dir) (fs/parent (first vcfs)) out-dir)
          file-info {:out-vcf (str (fs/file base-dir
                                            (itx/add-file-part (-> vcfs first fs/base-name)
                                                               (case merge-type
                                                                     :minimal "mincombine"
                                                                     :full "fullcombine"
                                                                     "combine"))))}
          args (concat ["-R" ref
                        "-o" :out-vcf
                        "--rod_priority_list" (join "," (map unique-name vcfs))]
                       (flatten (map #(list (str "--variant:" (unique-name %)) %) vcfs))
                       (if-not (nil? intervals) ["-L", intervals] [])
                       (case merge-type
                             :full ["--genotypemergeoption" "PRIORITIZE"]
                             :unique ["--genotypemergeoption" "UNIQUIFY"]
                             :minimal ["--sites_only" "--minimalVCF"]))]
      (if-not (fs/exists? base-dir)
        (fs/mkdirs base-dir))
      (broad/run-gatk "CombineVariants" args file-info {:out [:out-vcf]})
      (:out-vcf file-info))))

(defn convert-no-calls [in-vcf align-bam ref & {:keys [out-dir] :or {out-dir nil}}]
  "Convert no-calls into callable reference and real no-calls."
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

(defn select-by-sample [sample call ref & {:keys [out-dir intervals]
                                           :or {out-dir nil
                                                intervals nil}}]
  "Select only the sample of interest from input VCF files."
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

(defn create-merged [vcfs align-bams do-merges ref & {:keys [out-dir intervals]
                                                      :or {out-dir nil
                                                           intervals nil}}]
  "Create merged VCF files with no-call/ref-calls for each of the inputs."
  (letfn [(merge-vcf [vcf all-vcf align-bam ref]
            (let [ready-vcf (combine-variants [vcf all-vcf] ref
                                              :merge-type :full)]
              (convert-no-calls ready-vcf align-bam ref :out-dir out-dir)))]
    (let [merged (combine-variants vcfs ref :merge-type :minimal :intervals intervals)]
      (map (fn [[v b merge?]] (if merge? (merge-vcf v merged b ref) v))
           (map vector vcfs align-bams do-merges)))))

(defn gatk-normalize [call ref out-dir]
  "Prepare call information for VCF comparisons by normalizing through GATK.
  Handles:
   1. Combining multiple input files
   2. Splitting combined MNPs into phased SNPs"
  (if-not (fs/exists? out-dir)
    (fs/mkdirs out-dir))
  (letfn [(merge-call-files [call]
            (combine-variants (:file call) ref
                              :merge-type :full :out-dir out-dir))]
    (let [merge-file (if (coll? (:file call))
                       (merge-call-files call)
                       (:file call))]
      (assoc call :file (normalize-variants merge-file ref out-dir)))))
