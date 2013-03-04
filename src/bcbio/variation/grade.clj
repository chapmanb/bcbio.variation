(ns bcbio.variation.grade
  "Support comparisons of variant calls to reference call sets, providing
   detailed metrics about problematic discordant calls."
  (:import [org.broadinstitute.sting.utils.codecs.vcf
            VCFInfoHeaderLine VCFHeaderLineType]
           [org.broadinstitute.sting.utils.variantcontext VariantContextBuilder])
  (:require [clojure.set :refer [intersection]]
            [lonocloud.synthread :as ->]
            [bcbio.run.itx :as itx]
            [bcbio.variation.combine :as combine]
            [bcbio.variation.filter.attr :as attr]
            [bcbio.variation.variantcontext :as gvc]))

(defn- find-grading-and-eval
  "Separate grading reference and evaluation genome based on `type` parameter."
  [& cs]
  (let [type-groups (group-by #(keyword (:type %)) cs)]
    [(first (:grading-ref type-groups))
     (first (get type-groups nil))]))

(defn- to-refcalls
  "Convert truth discordants into reference calls "
  [f ref-file]
  (let [out-file (itx/add-file-part f "asref")]
    (when (itx/needs-run? out-file)
      (with-open [in-vcf-iter (gvc/get-vcf-iterator f ref-file)]
        (gvc/write-vcf-w-template f {:out out-file}
                                  (map #(gvc/genotypes->refcall % :attrs #{"DP"})
                                       (gvc/parse-vcf in-vcf-iter))
                                  ref-file)))
    out-file))

(defn- merge-discordants
  "Merge extra and missing discordant calls into single VCF."
  [eval-vcf truth-vcf ref-file]
  (combine/combine-variants [eval-vcf (to-refcalls truth-vcf ref-file)]
                            ref-file :merge-type :full :quiet-out? true))

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
  (let [pls (when ref-vc
              (attr/get-pls ref-vc))]
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
  (let [[truth-c eval-c] (find-grading-and-eval (:c1 cmp) (:c2 cmp))
        eval-kw (keyword (str (:name eval-c) "-discordant"))
        truth-kw (keyword (str (:name truth-c) "-discordant"))
        eval-vcf (get-in cmp [:c-files eval-kw])
        truth-vcf (get-in cmp [:c-files truth-kw])
        ref-file (get-in cmp [:exp :ref])
        base-eval-vcf (merge-discordants eval-vcf truth-vcf ref-file)
        out-vcf (itx/add-file-part eval-vcf "annotate")]
    (when (itx/needs-run? out-vcf)
      (with-open [ref-get (gvc/get-vcf-retriever ref-file truth-vcf)
                  eval-iter (gvc/get-vcf-iterator base-eval-vcf ref-file)]
        (gvc/write-vcf-w-template base-eval-vcf {:out out-vcf}
                                  (map (partial add-grade-cat ref-get)
                                       (gvc/parse-vcf eval-iter))
                                  ref-file :header-update-fn add-grade-header)))
    (assoc-in cmp [:c-files eval-kw] out-vcf)))