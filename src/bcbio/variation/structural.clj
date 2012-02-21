(ns bcbio.variation.structural
  "Handle structural variations for larger insertions, deletions and
  genome rearrangements."
  (:use [bcbio.variation.variantcontext :only (parse-vcf)]))

(defn sv-type
  "Determine the type of a structural variant; nil if not structural variant.
  Expected types are:

    - DEL: Deletion
    - INS: Insertion
    - DUP: Duplication
    - INV: Inversion
    - BND: Breakpoint end; paired with second variant
    - CNV: Copy number variation"
  [vc]
  (let [min-indel 100]
    (letfn [(max-allele-size [vc]
              (apply max (map #(.length %) (cons (:ref-allele vc) (:alt-alleles vc)))))
            (indel-type [vc]
              (if (> (.length (:ref-allele vc)) 0) :DEL :INS))
            (alt-sv-type [vc]
              (->> (-> vc :alt-alleles first .getDisplayString)
                   (re-find #"^<(\w+)(:|>)" )
                   second
                   keyword))]
      (cond
       (and (= "INDEL" (:type vc))
            (> (max-allele-size vc) min-indel)) (indel-type vc)
       (= "SYMBOLIC" (:type vc)) (alt-sv-type vc)
       :else nil))))
