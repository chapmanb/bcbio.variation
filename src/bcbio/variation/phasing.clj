;; Support phased haplotype comparisons between variants

(ns bcbio.variation.phasing
  (:use [bcbio.variation.variantcontext :only [parse-vcf get-vcf-retriever]]))

(defn is-phased?
  "Check for phasing on a single genotype variant context."
  [vc]
  {:pre [(= 1 (count (:genotypes vc)))]}
  (-> vc :genotypes first :genotype .isPhased))

(defn parse-phased-haplotypes [in-vcf]
  "Separate phased haplotypes provided in diploid input genome."
  (lazy-seq
   (loop [vcs (parse-vcf in-vcf)
          cur-hap []
          all-haps []]
     ;; 3 conditions:
     ;; 1. Out of variants; add the current one to the list and return
     ;; 2. No current haplotype variants or phased with the previous variant:
     ;;    add to the current haplotype
     ;; 3. A new haplotype: add existing haplotype to list and create new
     (cond
      (nil? (first vcs)) (if (empty? cur-hap) all-haps (conj all-haps cur-hap))
      (or (empty? cur-hap)
          (is-phased? (first vcs))) (recur (rest vcs) (conj cur-hap (first vcs)) all-haps)
      :else (recur (rest vcs) [(first vcs)] (conj all-haps cur-hap))))))

(defn parse-haploid-reference [in-vcf]
  "Prepare haploid reference variants for comparisons to called variants."
  (let [vcf-fetch (get-vcf-retriever in-vcf)]
    (doseq [vc (vcf-fetch "chr22" 16 16)]
      (println vc))))
