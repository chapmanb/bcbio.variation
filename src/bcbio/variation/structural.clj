(ns bcbio.variation.structural
  "Handle structural variations for larger insertions, deletions and
  genome rearrangements.")

(defn sv-type
  "Determine the type of a structural variant. Expected types are:

    - DEL: Deletion
    - INS: Insertion
    - DUP: Duplication
    - INV: Inversion
    - BND: Breakpoint end; paired with second variant
    - CNV: Copy number variation
    - nil: Not a structural variant."
  [vc]
  (let [min-indel 100]
    (letfn [(max-allele-size [vc]
              (apply max (map #(.length %) (cons (:ref-allele vc) (:alt-alleles vc)))))
            (indel-type [vc]
              (if (> (.length (:ref-allele vc)) 0) :DEL :INS))
            (sv-type-from-symbol [allele]
              (->> allele
                   (re-find #"^<(\w+)(:|>)" )
                   second
                   keyword))
            (alt-sv-type [vc]
              (let [allele (-> vc :alt-alleles first .getDisplayString)]
                (cond
                 (.startsWith allele "<") (sv-type-from-symbol allele)
                 (or (.contains allele "[")
                     (.contains allele "]")) :BND)))]
      (cond
       (and (= "INDEL" (:type vc))
            (> (max-allele-size vc) min-indel)) (indel-type vc)
       (= "SYMBOLIC" (:type vc)) (alt-sv-type vc)
       :else nil))))
