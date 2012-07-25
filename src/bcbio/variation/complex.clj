(ns bcbio.variation.complex
  "Handle complex variations representations: multi-nucleotide
   polymorphisms and indels."
  (:import [org.broadinstitute.sting.utils.variantcontext Allele
            VariantContextBuilder GenotypesContext GenotypeBuilder
            VariantContextUtils]
           [org.biojava3.core.sequence DNASequence]
           [org.biojava3.alignment Alignments])
  (:use [clojure.java.io]
        [clojure.set :only [union]]
        [ordered.set :only [ordered-set]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template
                                               get-vcf-source]])
  (:require [clojure.string :as string]
            [bcbio.run.broad :as broad]
            [bcbio.run.itx :as itx]))

;; ## Multi-nucleotide polymorphisms (MNPs)
;; Split into single variant primitives.

(defn- has-variant-base?
  "Do a set of alleles have any variants at a position."
  [alleles i]
  (> (count (set (map #(nth % i nil) alleles)))
     1))

(defn- get-vc-alleles [vc]
  (map #(.getBaseString %) (cons (:ref-allele vc) (:alt-alleles vc))))

(defn is-multi-indel?
  "Identify complex indels that can be split into multiple calls."
  [vc]
  (letfn [(monomorphic-alleles? [vc]
            (= 1 (->> (get-vc-alleles vc)
                      (map set)
                      (apply union)
                      count)))]
    (and (= "INDEL" (:type vc))
         (and (> (.length (:ref-allele vc)) 1)
              (> (apply max (map #(.length %) (:alt-alleles vc))) 1)
              (not (monomorphic-alleles? vc))))))

(defn- split-alleles
  "Detect single call SNP variants within a MNP genotype.
  Handles reference no-variant padding bases on the 5' end of
  the sequence, writing only variants at the adjusted positions."
  [vc alleles]
  (letfn [(mnp-ref-padding [ref-allele vc]
            {:post [(>= % 0)]}
            (- (inc (- (:end vc) (:start vc)))
               (count (string/replace ref-allele "-" ""))))
          (contains-indel? [alleles i]
            (contains? (set (map #(str (nth % i)) alleles)) "-"))
          (is-start-indel? [alleles start pos]
            (and (= 0 start)
                 (some empty? (map #(string/replace (subs % start pos) "-" "") alleles))))
          (extend-indels [alleles i]
            (if-not (contains-indel? alleles i) 
              {:start i :end (inc i)}
              {:start (max (dec i) 0)
               :end (inc (last (take-while #(or (contains-indel? alleles %) (is-start-indel? alleles i %))
                                           (range i (count (first alleles))))))}))
          (extract-variants [alleles pos ref-pad]
            (let [{:keys [start end]} (extend-indels alleles pos)
                  str-alleles (map #(-> (str (subs % start end))
                                        (string/replace "-" ""))
                                   alleles)
                  cur-alleles (map-indexed (fn [i x] (Allele/create x (= 0 i)))
                                           (into (ordered-set) str-alleles))]
              {:offset (+ start ref-pad)
               :end end
               :size (dec (.length (first cur-alleles)))
               :ref-allele (first cur-alleles)
               :alleles (rest cur-alleles)}))]
    (let [ref-pad (mnp-ref-padding (first alleles) vc)]
      (remove nil?
              (loop [i 0 final []]
                (cond
                 (> i (-> alleles first count)) final
                 (has-variant-base? alleles i)
                 (let [next-var (extract-variants alleles i ref-pad)]
                   (recur (:end next-var) (conj final next-var)))
                 :else (recur (inc i) final)))))))

(defn genotype-w-alleles
  "Retrieve a new genotype with the given alleles.
   Creates a single genotype from the VariantContext, copying the existing
   genotype and substituting in the provided alleles and phasing information."
  [vc alleles is-phased]
  (let [genotype (first (.getGenotypes vc))]
    (doto (-> vc .getGenotypes GenotypesContext/copy)
      (.replace
       (-> (GenotypeBuilder. genotype)
           (.alleles alleles)
           (.phased is-phased)
           .make)))))

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
  {:pre [(= 1 (:num-samples vc))]}
  (let [alleles (split-alleles vc (get-vc-alleles vc))]
    (map (fn [[i x]] (new-split-vc (:vc vc) i x)) (map-indexed vector alleles))))

;; ## Indels
;; Create a normalized representation for comparison.

(defn- multiple-alignment
  "Perform alignment of input sequences using BioJava."
  [seqs]
  (map #(.getSequenceAsString %)
       (-> (map #(DNASequence. %) seqs)
           (Alignments/getMultipleSequenceAlignment (to-array []))
           .getAlignedSequences)))

(defn- fix-ref-alignment-gaps
  "Ensure reference alignment gaps have consistent gap schemes."
  [alleles]
  (letfn [(has-5-gaps? [x] (.startsWith x "-"))
          (has-3-gaps? [x] (.endsWith x "-"))
          (has-internal-gaps? [x] (and (> (count x) 2)
                                       (.contains (subs x 1 (dec (count x))) "-")))
          (make-3-gap-only [x]
            (let [nogap-x (string/replace x "-" "")]
              (string/join "" (cons nogap-x
                                    (repeat (- (count x) (count nogap-x)) "-")))))]
    (let [ref-allele (first alleles)]
      (cons (if (or (has-internal-gaps? ref-allele)
                    (and (has-5-gaps? ref-allele) (has-3-gaps? (second alleles))))
              (make-3-gap-only ref-allele)
              ref-allele)
            (map #(if (has-internal-gaps? %) (make-3-gap-only %) %) (rest alleles))))))

(defn- split-complex-indel
  "Split complex indels into individual variant components."
  [vc]
  {:pre [(= 1 (:num-samples vc))]}
  (let [alleles (split-alleles vc (->> (get-vc-alleles vc)
                                       (remove empty?)
                                       multiple-alignment
                                       fix-ref-alignment-gaps))]
    (map (fn [[i x]] (new-split-vc (:vc vc) i x)) (map-indexed vector alleles))))

(defn- maybe-strip-indel
  "Remove extra variant bases, if necessary, from 5' end of indels.
  Checks both called alleles and potential alleles for extra 5' padding
  removing this if not needed to distinguish any potential alleles."
  [vc]
  {:pre [(= 1 (:num-samples vc))]}
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
  [vc-iter mnp-blockers]
  (letfn [(overlaps-previous-mnp? [vc blockers]
            (contains? (get blockers (:chr vc) #{}) (:start vc)))
          (add-mnp-info [coll vc]
            (let [prev-check 10000] ;; bp to check upstream for overlaps
              (if (or (= "MNP" (:type vc)) (is-multi-indel? vc))
                (assoc coll (:chr vc)
                       (union (set (remove #(< % (- (:start vc) prev-check))
                                           (get coll (:chr vc) #{})))
                              (set (range (:start vc)
                                          (+ (:start vc) (apply max (map count (get-vc-alleles vc))))))))
                coll)))
          (process-vc [vc blockers]
            (if (overlaps-previous-mnp? vc blockers)
              []
              (condp = (:type vc)
                "MNP" (split-mnp vc)
                "INDEL" (if (is-multi-indel? vc)
                          (split-complex-indel vc)
                          [(maybe-strip-indel vc)])
                [(:vc vc)])))
          (add-normalized-vcs [vcs blockers]
            (when (seq vcs)
              (concat (process-vc (first vcs) blockers)
                      (get-normalized-vcs (rest vcs) (add-mnp-info blockers (first vcs))))))]
    (lazy-seq (add-normalized-vcs vc-iter mnp-blockers))))

(defn left-align-variants
  "Left align variants in an input VCF file for a standard representation.
  Checks final line count of prepared file, returning left-aligned file
  when converting every variant in the input."
  [in-file ref & {:keys [out-dir]}]
  (letfn [(line-count [f]
            (with-open [rdr (reader f)]
              (count (line-seq rdr))))]
    (let [file-info {:out-vcf (itx/add-file-part in-file "leftalign" out-dir)}
          args ["-R" ref "-o" :out-vcf "--variant" in-file]]
      (broad/run-gatk "LeftAlignVariants" args file-info {:out [:out-vcf]})
      (if (= (line-count in-file) (line-count (:out-vcf file-info)))
        (:out-vcf file-info)
        in-file))))

(defn normalize-variants
  "Convert MNPs and indels into normalized representation."
  ([in-file ref]
     (normalize-variants in-file ref nil))
  ([in-file ref out-dir & {:keys [out-fname]}]
     (let [base-name (if (nil? out-fname) (itx/remove-zip-ext in-file) out-fname)
           out-file (itx/add-file-part base-name "nomnp" out-dir)]
       (when (itx/needs-run? out-file)
         (let [la-file (left-align-variants in-file ref :out-dir out-dir)]
           (with-open [vcf-source (get-vcf-source la-file ref)]
             (write-vcf-w-template in-file {:out out-file}
                                   (get-normalized-vcs (parse-vcf vcf-source) {})
                                   ref))))
       out-file)))
