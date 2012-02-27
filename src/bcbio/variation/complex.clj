(ns bcbio.variation.complex
  "Handle complex variations representations: multi-nucleotide
   polymorphisms and indels."
  (:import [org.broadinstitute.sting.utils.variantcontext Allele
            VariantContextBuilder GenotypesContext Genotype
            VariantContextUtils])
  (:use [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template]])
  (:require [bcbio.run.itx :as itx]))

;; ## Multi-nucleotide polymorphisms (MNPs)
;; Split into single variant primitives.

(defn- has-variant-base?
  "Do a set of alleles have any variants at a position."
  [alleles i]
  (> (count (set (map #(nth % i nil) alleles)))
     1))

(defn- split-alleles
  "Detect single call SNP variants within a MNP genotype."
  [vc genotype]
  (letfn [(extract-variants [alleles i]
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

(defn genotype-w-alleles
  "Retrieve a new genotype with the given alleles.
   Creates a single genotype from the VariantContext, copying the existing
   genotype and substituting in the provided alleles and phasing information."
  [vc alleles is-phased]
  (let [genotype (first (.getGenotypes vc))]
    (doto (-> vc .getGenotypes GenotypesContext/copy)
      (.replace
       (Genotype. (.getSampleName genotype)
                  alleles
                  (.getLog10PError genotype)
                  (if (.filtersWereApplied genotype) (.getFilters genotype) nil)
                  (.getAttributes genotype)
                  is-phased)))))

(defn- new-split-vc
  "Create a new VariantContext as a subset of an existing variant.
   `allele-info` specifies the location size and alleles for the new variant:
   `{:offset :size :ref-allele :alleles}`"
  [vc i allele-info]
  (let [pos (+ (:offset allele-info) (.getStart vc))]
    (-> (VariantContextBuilder. vc)
        (.start pos)
        (.stop (+ pos (get allele-info :size 0)))
        (.genotypes (genotype-w-alleles vc (:alleles allele-info) (> i 0)))
        (.alleles (set (cons (:ref-allele allele-info) (:alleles allele-info))))
        (.make))))

(defn- split-mnp
  "Split a MNP into individual alleles"
  [vc]
  {:pre [(= 1 (count (:genotypes vc)))]}
  (let [alleles (split-alleles vc (-> vc :genotypes first))]
    (map (fn [[i x]] (new-split-vc (:vc vc) i x)) (map-indexed vector alleles))))

;; ## Indels
;; Create a normalized representation for comparison.

(defn- maybe-strip-indel
  "Remove extra variant bases, if necessary, from 5' end of indels."
  [vc]
  {:pre [(= 1 (count (:genotypes vc)))]}
  (letfn [(strip-indel [vc i alleles]
            (let [start-pos (- i 1)
                  ref-allele (subs (first alleles) start-pos)
                  cur-alleles (map #(Allele/create (subs % start-pos)
                                                   (= ref-allele (subs % start-pos)))
                                   alleles)]
              (new-split-vc vc 0 {:offset i
                                  :size (- (count ref-allele) 1)
                                  :ref-allele (first cur-alleles)
                                  :alleles (rest cur-alleles)})))]
    (let [orig-alleles (map #(.getBaseString %) (cons (:ref-allele vc)
                                                      (-> vc :genotypes first :alleles)))
          first-var-i (first (filter #(has-variant-base? orig-alleles %)
                                     (range (apply max (map count orig-alleles)))))]
      (if (or (nil? first-var-i) (= first-var-i 0))
        (:vc vc)
        (strip-indel (:vc vc) first-var-i orig-alleles)))))

;; ## VCF file conversion
;; Process entire files, normalizing complex variations

(defn- get-normalized-vcs
  "Lazy list of variant context with MNPs split into single genotypes and indels stripped."
  [in-file]
  (map (fn [x] [:out x])
   (flatten
    (for [vc (parse-vcf in-file)]
      (condp = (:type vc)
        "MNP" (split-mnp vc)
        "INDEL" (maybe-strip-indel vc)
        (:vc vc))))))

(defn normalize-variants
  "Convert MNPs and indels into normalized representation."
  ([in-file ref]
     (normalize-variants in-file ref nil))
  ([in-file ref out-dir & {:keys [out-fname]}]
     (let [base-name (if (nil? out-fname) (itx/remove-zip-ext in-file) out-fname)
           out-file (itx/add-file-part base-name "nomnp" out-dir)]
       (if (itx/needs-run? out-file)
         (write-vcf-w-template in-file {:out out-file} (get-normalized-vcs in-file) ref))
       out-file)))
