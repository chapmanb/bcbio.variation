;; Support phased haplotype comparisons between variants

(ns bcbio.variation.phasing
  (:use [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template]]))

(defn is-phased?
  "Check for phasing on a single genotype variant context."
  [vc]
  {:pre [(= 1 (count (:genotypes vc)))]}
  (-> vc :genotypes first :genotype .isPhased))

(defn parse-phased-haplotypes [in-vcf]
  "Separate phased haplotypes provided in diploid input genome."
  (partition-by is-phased? (parse-vcf in-vcf)))
