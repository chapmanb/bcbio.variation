(ns bcbio.variation.grade
  "Support comparisons of variant calls to reference call sets, providing
   detailed metrics about problematic discordant calls."
  (:import [org.broadinstitute.variant.vcf
            VCFInfoHeaderLine VCFHeaderLineType]
           [org.broadinstitute.variant.variantcontext VariantContextBuilder])
  (:require [clojure.set :refer [intersection]]
            [clojure.string :as string]
            [clojure.math.combinatorics :as combo]
            [lonocloud.synthread :as ->]
            [bcbio.run.itx :as itx]
            [bcbio.variation.annotation :as annotation]
            [bcbio.variation.combine :as combine]
            [bcbio.variation.filter.attr :as attr]
            [bcbio.variation.phasing :as phasing]
            [bcbio.variation.variantcontext :as gvc]))

;; ## Utility functions

(defn is-grade-cmp?
  [exp]
  (= :grade (keyword (get exp :approach "compare"))))

(defn is-grading-ref?
  "Grading references are either haploid or identified with grading-ref type."
  [exp c]
  (or (= :grading-ref (keyword (:type c)))
      (-> c :file (phasing/is-haploid? (:ref exp)))))

(defn- find-grading-and-eval-kws
  "Separate grading reference and evaluation genome based on `type` parameter.
   Defaults to the first being reference and second evaluation if not defined."
  [exp c1 c2 c-files]
  (let [grade-groups (group-by (partial is-grading-ref? exp) [c1 c2])
        marked-ref (first (get grade-groups true))
        marked-eval (first (get grade-groups false))
        [truth-c eval-c] (if (and marked-ref marked-eval)
                           [marked-ref marked-eval]
                           [c1 c2])
        eval-kw (keyword (str (:name eval-c) "-discordant"))
        truth-kw (keyword (str (:name truth-c) "-discordant"))]
    {:eval (if (contains? c-files eval-kw) eval-kw :discordant)
     :truth (if (contains? c-files truth-kw) truth-kw :discordant-missing)}))

;; ## Summarize grading results

(defn- pick-discordant-reason
  "Decide on a likely reason for a discordant variant call"
  [vc attr-getter]
  (letfn [(is-repeat-region? [attrs]
            (or (< (or (get attrs "gms_illumina") 100.0) 50.0)
                (contains? (or (get attrs "rmsk") #{}) "repeat")))
          (is-error-prone? [attrs]
            (contains? (or (get attrs "in_cse") #{}) "error-prone"))]
    (let [attrs (attr-getter ["DP" "rmsk" "gms_illumina" "in_cse"] vc)]
      (cond
       (< (or (get attrs "DP") 500) 10) :low-coverage
       (is-error-prone? attrs) :error-prone
       (is-repeat-region? attrs) :repeat
       :else :other))))

(defn- identify-discordant-cat
  "Identify the variant type and discordant category.
   - variant types -- :snp :indel
   - discordant types -- :shared :missing :extra
   - reason types -- :hethom :vardiff :low-coverage :repeat :error-prone :other"
  [vc attr-getter]
  (let [vtype (keyword (string/lower-case (:type vc)))
        cat (-> (get-in vc [:attributes "GradeCat"])
                (string/replace "discordant-" "")
                keyword)
        [dtype rtype] (case cat
                        (:missing :extra) [cat (pick-discordant-reason vc attr-getter)]
                        [:shared cat])]
    [vtype dtype rtype]))

(defn- count-discordant-categories
  [vcf-file ref-file]
  (let [attr-getter (attr/prep-vc-attr-retriever vcf-file ref-file)]
    (with-open [in-vcf-iter (gvc/get-vcf-iterator vcf-file ref-file)]
      (reduce (fn [coll vc]
                (let [cat (identify-discordant-cat vc attr-getter)]
                  (assoc-in coll cat (inc (get-in coll cat 0)))))
              {} (gvc/parse-vcf in-vcf-iter)))))

(defn prep-grade-breakdown
  "Prepare detailed grading breakdown of concordant and discordant variants.
   The goal is to help identify common causes of discordance."
  [cmp]
  (let [kws (find-grading-and-eval-kws (:exp cmp) (:c1 cmp) (:c2 cmp)
                                       (:c-files cmp))
        vcf-file (get-in cmp [:c-files (:eval kws)])
        ref-file (get-in cmp [:exp :ref])
        summary (:summary cmp)]
    {:sample (:sample summary)
     :concordant (select-keys summary [:genotype_concordance :callable_concordance
                                       :concordant])
     :discordant (count-discordant-categories vcf-file ref-file)}))

;; ## Identify grading references

(defn- to-refcalls
  "Convert truth discordants into reference calls "
  [f ref-file]
  (let [out-file (itx/add-file-part f "asref")]
    (when (itx/needs-run? out-file)
      (with-open [in-vcf-iter (gvc/get-vcf-iterator f ref-file)]
        (gvc/write-vcf-w-template f {:out out-file}
                                  (map #(gvc/genotypes->refcall %)
                                       (gvc/parse-vcf in-vcf-iter))
                                  ref-file)))
    out-file))

(defn- merge-discordants
  "Merge extra and missing discordant calls into single VCF."
  [eval-vcf truth-vcf align-bam ref-file]
  (let [truth-dis-vcf (-> truth-vcf
                          (to-refcalls ref-file)
                          (->/when align-bam
                            (annotation/add-gatk-annotations align-bam ref-file :annos ["DepthOfCoverage"])))]
    (combine/combine-variants [eval-vcf truth-dis-vcf]
                              ref-file :merge-type :full :quiet-out? true :check-ploidy? false)))

;; ## Add grading info to VCF

(defn hethom-discordant?
  "Identify variants that are variant calls but discordant based on het/hom calls."
  [vc other-vcs]
  (letfn [(get-alleles [x]
            (->> (:genotypes x)
                 (map :alleles)
                 flatten
                 (remove #(.isReference %))
                 set))]
    (let [vc2 (first (filter #(and (= (:start vc) (:start %))
                                   (= (:ref-allele vc) (:ref-allele %)))
                             other-vcs))]
      (seq (intersection (get-alleles vc) (get-alleles vc2))))))

(defn is-novar-call?
  "Is a variant context a no-variant (reference) call."
  [vc]
  (->> (:genotypes vc)
       (map :alleles)
       flatten
       (every? #(.isReference %))))

(defn- assign-grade-cat
  "Assign a discordance category to variants that do not match grading reference."
  [vc ref-vcs]
  (cond
   (empty? ref-vcs) :discordant-extra
   (is-novar-call? vc) :discordant-missing
   (hethom-discordant? vc ref-vcs) :discordant-hethom
   :else :discordant-vardiff))

(defn- get-grade-score
  "Determine likelihood of grading reference based on grading category.
   :discordant-missing -- probability that grading standard is actually reference.
   :discordant-hethom  -- probability that grading standard is alternative variant."
  [ref-vc cat]
  {:pre [(or (nil? ref-vc) (= 1 (:num-samples ref-vc)))]}
  (let [pls (when ref-vc
              (attr/get-pls (-> ref-vc :genotypes first)))]
    (case cat
      :discordant-missing (get pls "HOM_REF")
      :discordant-hethom (first (vals (dissoc pls "HOM_REF")))
      nil)))

(defn add-grade-cat
  "Add grading category and score, providing additional details on discordant variants."
  [ref-get vc]
  (let [ref-vcs (gvc/variants-in-region ref-get vc)
        grade-cat (assign-grade-cat vc ref-vcs)]
    (-> (VariantContextBuilder. (:vc vc))
        (.attributes (-> (:attributes vc)
                         (assoc "GradeCat" (name grade-cat))
                         (->/when-let [score (get-grade-score (first ref-vcs) grade-cat)]
                           (assoc "GradeScore" (float score)))))
        .make)))

(defn add-grade-header
  "Add grading INFO fields to VCF output header"
  [_ header]
  (gvc/header-w-md
   header
   #{(VCFInfoHeaderLine. "GradeCat" 1 VCFHeaderLineType/String
                         "Grading category based on comparison with reference call set.")
     (VCFInfoHeaderLine. "GradeScore" 1 VCFHeaderLineType/Float
                         "Grading score: phred score of correct probability reference call.")}))

(defn annotate-discordant
  "Update a comparison with annotated information on discordant grading"
  [cmp]
  (let [kws (find-grading-and-eval-kws (:exp cmp) (:c1 cmp) (:c2 cmp)
                                       (:c-files cmp))
        eval-vcf (get-in cmp [:c-files (:eval kws)])
        truth-vcf (get-in cmp [:c-files (:truth kws)])
        ref-file (get-in cmp [:exp :ref])
        align-bam (get-in cmp [:exp :align])
        base-eval-vcf (merge-discordants eval-vcf truth-vcf align-bam ref-file)
        out-vcf (itx/add-file-part eval-vcf "annotate")]
    (when (itx/needs-run? out-vcf)
      (with-open [ref-get (gvc/get-vcf-retriever ref-file truth-vcf)
                  eval-iter (gvc/get-vcf-iterator base-eval-vcf ref-file)]
        (gvc/write-vcf-w-template base-eval-vcf {:out out-vcf}
                                  (map (partial add-grade-cat ref-get)
                                       (gvc/parse-vcf eval-iter))
                                  ref-file :header-update-fn add-grade-header)))
    (assoc-in cmp [:c-files (:eval kws)] out-vcf)))
