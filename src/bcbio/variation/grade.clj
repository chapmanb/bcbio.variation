(ns bcbio.variation.grade
  "Support comparisons of variant calls to reference call sets, providing
   detailed metrics about problematic discordant calls."
  (:require [bcbio.run.itx :as itx]
            [bcbio.variation.combine :as combine]
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
                  eval-iter (gvc/get-vcf-iterator base-eval-vcf ref-file)]))
    (assoc-in cmp [:c-files eval-kw] out-vcf)))