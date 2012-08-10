(ns bcbio.variation.complex
  "Handle complex variations representations: multi-nucleotide
   polymorphisms and indels."
  (:import [org.broadinstitute.sting.utils.variantcontext Allele
            VariantContextBuilder GenotypesContext GenotypeBuilder
            VariantContextUtils]
           [org.biojava3.core.sequence DNASequence]
           [org.biojava3.alignment Alignments SimpleGapPenalty
            Alignments$PairwiseSequenceScorerType])
  (:use [clojure.java.io]
        [clojure.set :only [union]]
        [ordered.set :only [ordered-set]]
        [bcbio.align.ref :only [extract-sequence]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template
                                               get-vcf-iterator pad-vc-alleles]])
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
  (vec (map #(.getBaseString %) (cons (:ref-allele vc) (:alt-alleles vc)))))

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
              (> (apply min (map #(.length %) (:alt-alleles vc))) 1)
              (not (monomorphic-alleles? vc))))))

(defn- split-alleles
  "Detect single call SNP variants within a MNP genotype.
  Handles reference no-variant padding bases on the 5' end of
  the sequence, writing only variants at the adjusted positions."
  [vc alleles]
  (letfn [(contains-indel? [alleles i]
            (when (< i (count (first alleles)))
              (contains? (set (map #(str (nth % i)) alleles)) "-")))
          (starts-an-indel? [alleles i]
            (contains-indel? alleles (inc i)))
          (is-match? [alleles i]
            (= 1 (count (set (map #(str (nth % i)) alleles)))))
          (is-internal-indel? [alleles i]
            (and (pos? i)
                 (contains-indel? alleles i)
                 (is-match? alleles (dec i))))
          (is-fiveprime-indel? [alleles i]
            (and (zero? i)
                 (or
                  (starts-an-indel? alleles i)
                  (contains-indel? alleles i))))
          (needs-padding? [alleles i]
            (or (pos? i)
                (and (is-fiveprime-indel? alleles i)
                     (is-match? alleles i))
                (and (not (is-fiveprime-indel? alleles i))
                     (not (is-match? alleles i)))))
          (mnp-ref-padding [alleles vc]
            {:post [(>= % 0)]}
            (let [first-no-indel (first (drop-while #(contains-indel? alleles %)
                                                    (range 0 (count (first alleles)))))
                  gatk-pad (max 0 (- (inc (- (:end vc) (:start vc)))
                                     (count (string/replace (first alleles) "-" ""))))]
              (cond (or (is-match? alleles 0) (pos? gatk-pad)) gatk-pad
                    (and (contains-indel? alleles 0)
                         (not (is-match? alleles first-no-indel))) 0
                    :else 0)))
          (extend-indels [alleles i]
            {:start (if (or (is-internal-indel? alleles i)
                            (is-fiveprime-indel? alleles i))
                      (max (dec i) 0)
                      i)
             :end (inc (or (last (take-while #(or (contains-indel? alleles %)
                                                  (starts-an-indel? alleles %))
                                             (range i (count (first alleles)))))
                           i))})
          (extract-variants [alleles pos ref-pad]
            (let [{:keys [start end]} (extend-indels alleles pos)
                  str-alleles (map #(-> (str (subs % start end))
                                        (string/replace "-" ""))
                                   alleles)
                  cur-alleles (map-indexed (fn [i x] (Allele/create x (= 0 i)))
                                           (into (ordered-set) str-alleles))
                  size (let [base (.length (first cur-alleles))]
                         (if (some empty? str-alleles) base (dec base)))
                  w-gap-start (-> (first alleles) (subs 0 start) (string/replace "-" "") count)]
              {:offset (+ w-gap-start 
                          (cond (needs-padding? alleles start) ref-pad
                                (and (is-fiveprime-indel? alleles start)
                                     (contains-indel? alleles start)) -1
                                :else 0))
               :end (+ ref-pad w-gap-start size)
               :next-start end
               :size size
               :ref-allele (first cur-alleles)
               :alleles (rest cur-alleles)}))]
    (let [ref-pad (mnp-ref-padding alleles vc)]
      (remove nil?
              (loop [i 0 final []]
                (cond
                 (>= i (-> alleles first count)) final
                 (has-variant-base? alleles i)
                 (let [next-var (extract-variants alleles i ref-pad)]
                   (recur (:next-start next-var) (conj final next-var)))
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
  (letfn [(original-seq-position [seqs]
            (let [orig-order (into {} (map-indexed (fn [i x] [x i])
                                                   (into (ordered-set) seqs)))]
              (fn [x]
                (get orig-order (string/replace x "-" "")))))
          (unique-aligns [xs]
            (vals (reduce (fn [coll x]
                            (assoc coll (string/replace x "-" "") x))
                          {} xs)))
          (all-gap? [xs]
            (= (set (map str xs)) #{"-"}))
          (finalize-alignment [seqs]
            (let [n 2
                  gap-free (remove all-gap? (partition n (apply interleave (take n seqs))))]
              [(string/join "" (map first gap-free))
               (string/join "" (map second gap-free))]))]
    (let [align-args (to-array [(SimpleGapPenalty. 20 1)])
          base-align (map #(.getSequenceAsString %)
                          (-> (map #(DNASequence. %) seqs)
                              (Alignments/getMultipleSequenceAlignment align-args)
                              .getAlignedSequences))
          orig-align (sort-by (original-seq-position seqs) (unique-aligns base-align))]
      (finalize-alignment orig-align))))

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
      (cons (if (has-internal-gaps? ref-allele) (make-3-gap-only ref-allele) ref-allele)
            (rest alleles)))))

(defn- sanity-check-split-vcs
  "Confirm that new variants match correctly back to original.
   Catch any potential errors in splitting by ensuring reference coordinates
   and sequences match original."
  [vc new-vcs]
  (letfn [(get-vc-info [vc]
            (let [alleles (map #(.getBaseString %) (.getAlleles vc))]
              {:start (.getStart vc)
               :pad (inc (- (- (.getEnd vc) (.getStart vc)) (count (first alleles))))
               :alleles alleles}))
          (check-split-vc [orig new]
            (let [int-pos (max 0 (- (+ (:start new) (:pad new))
                                    (+ (:start orig) (:pad orig))))
                  check-ref (first (:alleles new))]
              (when-not (or (>= int-pos (count (first (:alleles orig))))
                            (= (subs (first (:alleles orig)) int-pos (+ int-pos (count check-ref)))
                               check-ref))
                (throw (Exception. (format "Matching problem with split alleles: %s %s %s %s"
                                           (:chr vc) (:start vc) orig new))))))]
    (doall (map (partial check-split-vc (get-vc-info (:vc vc)))
                (map get-vc-info new-vcs)))))

(defn- split-complex-indel
  "Split complex indels into individual variant components."
  [orig-vc ref]
  {:pre [(= 1 (:num-samples orig-vc))]}
  (let [vc (pad-vc-alleles orig-vc)
        alleles (split-alleles vc (->> (conj (get-vc-alleles vc)
                                             (extract-sequence ref (:chr vc) (:start vc) (:end vc)))
                                       (remove empty?)
                                       (remove nil?)
                                       multiple-alignment))]
    (when-not (= (count alleles) (count (set (map :offset alleles))))
      (throw (Exception. (format "Mutiple alleles at same position: %s %s %s"
                                 (:chr vc) (:start vc) (vec alleles)))))
    (let [split-vcs (map (fn [[i x]] (new-split-vc (:vc vc) i x))
                         (map-indexed vector alleles))]
      (sanity-check-split-vcs vc split-vcs)
      split-vcs)))

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
  [vc-iter mnp-blockers ref]
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
                          (split-complex-indel vc ref)
                          [(maybe-strip-indel vc)])
                [(:vc vc)])))
          (add-normalized-vcs [vcs blockers]
            (when (seq vcs)
              (concat (process-vc (first vcs) blockers)
                      (get-normalized-vcs (rest vcs) (add-mnp-info blockers (first vcs)) ref))))]
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
           (with-open [vcf-iter (get-vcf-iterator la-file ref)]
             (write-vcf-w-template in-file {:out out-file}
                                   (get-normalized-vcs (parse-vcf vcf-iter) {} ref)
                                   ref))))
       out-file)))
