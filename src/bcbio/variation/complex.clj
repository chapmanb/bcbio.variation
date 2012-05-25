(ns bcbio.variation.complex
  "Handle complex variations representations: multi-nucleotide
   polymorphisms and indels."
  (:import [org.broadinstitute.sting.utils.variantcontext Allele
            VariantContextBuilder GenotypesContext Genotype
            VariantContextUtils])
  (:use [clojure.set :only [union]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template
                                               get-vcf-source]])
  (:require [bcbio.run.itx :as itx]))

;; ## Multi-nucleotide polymorphisms (MNPs)
;; Split into single variant primitives.

(defn- has-variant-base?
  "Do a set of alleles have any variants at a position."
  [alleles i]
  (> (count (set (map #(nth % i nil) alleles)))
     1))

(defn- split-alleles
  "Detect single call SNP variants within a MNP genotype.
  Handles reference no-variant padding bases on the 5' end of
  the sequence, writing only variants at the adjusted positions."
  [vc genotype]
  (letfn [(mnp-ref-padding [ref-allele vc]
            {:post [(>= % 0)]}
            (- (inc (- (:end vc) (:start vc)))
               (count ref-allele)))
          (extract-variants [alleles i ref-pad]
            (let [ref-allele (nth (first alleles) i)
                  cur-alleles (map #(Allele/create (str (nth % i))
                                                   (= ref-allele (nth % i)))
                                   alleles)]
              {:offset (+ i ref-pad)
               :ref-allele (first cur-alleles)
               :alleles (rest cur-alleles)}))]
    (let [orig-alleles (map #(.getBaseString %) (cons (:ref-allele vc) (:alleles genotype)))
          ref-pad (mnp-ref-padding (first orig-alleles) vc)]
      (remove nil?
              (for [i (-> orig-alleles first count range)]
                (if (has-variant-base? orig-alleles i)
                  (extract-variants orig-alleles i ref-pad)))))))

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
  "Remove extra variant bases, if necessary, from 5' end of indels.
  Checks both called alleles and potential alleles for extra 5' padding
  removing this if not needed to distinguish any potential alleles."
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
                                  :alleles (rest cur-alleles)})))
          (variant-allele-pos [input-alleles]
            (let [str-alleles (map #(.getBaseString %) input-alleles)
                  first-var-i (first (filter #(has-variant-base? str-alleles %)
                                     (range (apply max (map count str-alleles)))))]
              [str-alleles first-var-i]))]
    (let [[orig-alleles first-var-i] (variant-allele-pos (cons (:ref-allele vc)
                                                               (-> vc :genotypes first :alleles)))
          [_ nocall-i] (variant-allele-pos (cons (:ref-allele vc) (:alt-alleles vc)))]
      (if (or (nil? first-var-i) (= first-var-i 0)
              (nil? nocall-i) (= nocall-i 0))
        (:vc vc)
        (strip-indel (:vc vc) first-var-i orig-alleles)))))

;; ## VCF file conversion
;; Process entire files, normalizing complex variations

(defn- get-normalized-vcs
  "Lazy list of variant context with MNPs split into single genotypes and indels stripped.
  mnp-blockers is a lookup dictionary of positions overlapped by MNPs.
  Variant representations can have a MNP and also a single SNP representing
  the same information. In this case we ignore SNPs overlapping a MNP region
  and rely on MNP splitting to resolve the SNPs."
  [vcs mnp-blockers]
  (letfn [(overlaps-previous-mnp? [vc blockers]
            (contains? (get mnp-blockers (:chr vc) #{}) (:start vc)))
          (add-mnp-info [coll vc]
            (if (= "MNP" (:type vc))
              (assoc coll (:chr vc)
                     (union (get coll (:chr vc) #{})
                            (set (range (inc (:start vc)) (inc (:end vc))))))
              coll))
          (process-vc [vc blockers]
            (if (overlaps-previous-mnp? vc blockers)
              []
              (condp = (:type vc)
                "MNP" (split-mnp vc)
                "INDEL" [(maybe-strip-indel vc)]
                [(:vc vc)])))]
    (lazy-seq
     (when (seq vcs)
       (concat (process-vc (first vcs) mnp-blockers)
               (get-normalized-vcs (rest vcs) (add-mnp-info mnp-blockers (first vcs))))))))

(defn normalize-variants
  "Convert MNPs and indels into normalized representation."
  ([in-file ref]
     (normalize-variants in-file ref nil))
  ([in-file ref out-dir & {:keys [out-fname]}]
     (let [base-name (if (nil? out-fname) (itx/remove-zip-ext in-file) out-fname)
           out-file (itx/add-file-part base-name "nomnp" out-dir)]
       (when (itx/needs-run? out-file)
         (with-open [vcf-source (get-vcf-source in-file ref)]
           (write-vcf-w-template in-file {:out out-file}
                                 (get-normalized-vcs (parse-vcf vcf-source) {})
                                 ref)))
       out-file)))
