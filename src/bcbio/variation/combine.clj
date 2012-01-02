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
        [clojure.string :only [join]])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]))

(defn combine-variants [vcfs ref & {:keys [merge-type]
                                    :or {merge-type :unique}}]
  "Combine two variant files with GATK CombineVariants."
  (letfn [(unique-name [f]
            (-> f fs/base-name itx/file-root))]
    (let [file-info {:out-vcf (itx/add-file-part (first vcfs)
                                                 (case merge-type
                                                       :minimal "mincombine"
                                                       :full "fullcombine"
                                                       "combine"))}
          args (concat ["-R" ref
                        "-o" :out-vcf
                        "--rod_priority_list" (join "," (map unique-name vcfs))]
                       (flatten (map #(list (str "--variant:" (unique-name %)) %) vcfs))
                       (case merge-type
                             :full ["--genotypemergeoption" "REQUIRE_UNIQUE"]
                             :unique ["--genotypemergeoption" "UNIQUIFY"]
                             :minimal ["--sites_only" "--minimalVCF"]))]
      (broad/run-gatk "CombineVariants" args file-info {:out [:out-vcf]})
      (:out-vcf file-info))))

(defn convert-no-calls [in-vcf align-bam ref]
  "Convert no-calls into callable reference and real no-calls."
  (let [out-file (itx/add-file-part in-vcf "wrefs")
        is-callable? (callable-checker align-bam ref)]
    (letfn [(ref-genotype [g vc]
              (doto (-> vc :vc .getGenotypes GenotypesContext/copy)
                (.replace
                 (Genotype/modifyAlleles (:genotype g)
                                         (repeat (count (:alleles g))
                                                 (:ref-allele vc))))))
            (maybe-callable-vc [vc]
              {:pre (= 1 (count (:genotypes vc)))}
              (let [g (-> vc :genotypes first)]
                {:vc
                 (if (.isNoCall (-> g :alleles first))
                   (if (is-callable? (:chr vc) (:start vc) (:end vc))
                     (-> (VariantContextBuilder. (:vc vc))
                         (.genotypes (ref-genotype g vc))
                         (.make))
                     (-> (VariantContextBuilder. (:vc vc))
                         (.filters #{"NotCallable"})
                         (.make)))
                   (:vc vc))}))
            (convert-vcs [in-file]
              (for [vc (parse-vcf in-file)]
                [:out (maybe-callable-vc vc)]))]
      (if (itx/needs-run? out-file)
        (write-vcf-w-template in-vcf {:out out-file} (convert-vcs in-vcf) ref))
      out-file)))

(defn create-merged [vcfs align-bams ref]
  "Create merged VCF files with no-call/ref-calls for each of the inputs."
  (letfn [(merge-vcf [vcf all-vcf align-bam ref]
            (let [ready-vcf (combine-variants [vcf all-vcf] ref :merge-type :full)]
              (convert-no-calls ready-vcf align-bam ref)))]
    (let [merged (combine-variants vcfs ref :merge-type :minimal)]
      (map (fn [[v b]] (merge-vcf v merged b ref))
           (map vector vcfs align-bams)))))

