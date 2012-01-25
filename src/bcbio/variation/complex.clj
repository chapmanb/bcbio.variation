;; Handle more complex variations representations: multi-nucleotide
;; polymorphisms, indels and structural rearrangements

(ns bcbio.variation.complex
  (:import [org.broadinstitute.sting.utils.variantcontext Allele
            VariantContextBuilder GenotypesContext Genotype
            VariantContextUtils])
  (:use [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template]])
  (:require [bcbio.run.itx :as itx]))

;; Multi-nucleotide polymorphisms (MNPs): Split into single variant primitives.

(defn- split-alleles [vc genotype]
  "Split alleles for a MNP genotype into single call for each."
  (letfn [(has-variant-base? [alleles i]
            (> (count (set (map #(nth % i) alleles)))
               1))
          (extract-variants [alleles i]
            (let [ref-allele (nth (first alleles) i)
                  cur-alleles (map #(Allele/create (str (nth % i))
                                                   (= ref-allele (nth % i)))
                                   alleles)]
              {:offset i
               :ref-allele (first cur-alleles)
               :alleles (rest cur-alleles)}))]
    (let [orig-alleles (map #(.getBaseString %) (cons (:ref-allele vc) (:alleles genotype)))]
      (remove nil?
              (for [i (-> orig-alleles first count range)]
                (if (has-variant-base? orig-alleles i)
                  (extract-variants orig-alleles i)))))))

(defn- genotype-w-alleles [vc alleles is-phased]
  "Retrieve a new genotype with the given alleles."
  (let [genotype (first (.getGenotypes vc))]
    (doto (-> vc .getGenotypes GenotypesContext/copy)
      (.replace
       (Genotype. (.getSampleName genotype)
                  alleles
                  (.getLog10PError genotype)
                  (if (.filtersWereApplied genotype) (.getFilters genotype) nil)
                  (.getAttributes genotype)
                  is-phased)))))

(defn- new-split-vc [vc i allele-info]
  "Create a new VariantContext from a base MNP given a split allele."
  (let [pos (+ (:offset allele-info) (.getStart vc))]
    (-> (VariantContextBuilder. vc)
        (.start pos)
        (.stop pos)
        (.genotypes (genotype-w-alleles vc (:alleles allele-info) (> i 0)))
        (.alleles (set (cons (:ref-allele allele-info) (:alleles allele-info))))
        (.make))))

(defn- split-mnp
  "Split a MNP into individual alleles"
  [vc]
  {:pre [(= 1 (count (:genotypes vc)))]}
  (let [alleles (split-alleles vc (-> vc :genotypes first))]
    (map (fn [[i x]] (new-split-vc (:vc vc) i x)) (map-indexed vector alleles))))

;; Handle indels creating a normalized representation for comparison.

(defn- strip-indel [vc]
  "Remove extra variant bases "
  (-> vc
      ;(VariantContextUtils/createVariantContextWithPaddedAlleles true)
      (VariantContextUtils/createVariantContextWithTrimmedAlleles)))

(defn- get-normalized-vcs [in-file]
  "Lazy list of variant context with MNPs split into single genotypes and indels stripped."
  (map (fn [x] [:out x])
   (flatten
    (for [vc (parse-vcf in-file)]
      (condp = (:type vc)
        "MNP" (split-mnp vc)
        "INDEL" (strip-indel (:vc vc))
        (:vc vc))))))

(defn normalize-variants
  "Convert MNPs and indels into normalized representation."
  ([in-file ref]
     (normalize-variants in-file ref nil))
  ([in-file ref out-dir]
     (let [out-file (itx/add-file-part in-file "nomnp" out-dir)]
       (write-vcf-w-template in-file {:out out-file} (get-normalized-vcs in-file) ref)
       out-file)))
