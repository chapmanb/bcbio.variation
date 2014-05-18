(ns bcbio.variation.complex
  "Handle complex variations representations: multi-nucleotide
   polymorphisms and indels."
  (:import [org.broadinstitute.variant.variantcontext Allele
            VariantContextBuilder GenotypesContext GenotypeBuilder]
           [org.biojava3.core.sequence DNASequence]
           [org.biojava3.alignment Alignments SimpleGapPenalty
            Alignments$PairwiseSequenceScorerType])
  (:use [clojure.java.io]
        [clojure.set :only [union]]
        [ordered.set :only [ordered-set]]
        [bcbio.align.ref :only [extract-sequence]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template
                                               get-vcf-iterator]])
  (:require [clojure.string :as string]
            [bcbio.run.broad :as broad]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [me.raynes.fs :as fs]))

;; ## Multi-nucleotide polymorphisms (MNPs)
;; Split into single variant primitives.

(defn- has-variant-base?
  "Do a set of alleles have any variants at a position."
  [alleles i]
  (> (count (set (map #(nth % i nil) alleles)))
     1))

(defn- get-vc-alleles [vc]
  (vec (map #(.getDisplayString %) (cons (:ref-allele vc) (:alt-alleles vc)))))

(defn- is-multi-indel?
  "Identify complex indels that can be split into multiple calls.
   Caps indels we can operate on at 5000bp to avoid realignment errors for longer."
  [vc]
  (letfn [(monomorphic-alleles? [vc]
            (= 1 (->> (get-vc-alleles vc)
                      (map set)
                      (apply union)
                      count)))
          (has-multiple-nonref-alleles? [vc]
            (and (> (.length (:ref-allele vc)) 1)
                 (> (apply min (map #(.length %) (:alt-alleles vc))) 1)
                 (not (monomorphic-alleles? vc))))
          (has-ref-padding-mismatch? [vc]
            (let [alleles (get-vc-alleles vc)]
              (not= (nth (first alleles) 0) (nth (second alleles) 0))))
          (splittable-size? [vc]
            (< (apply max (map count (get-vc-alleles vc)))
               5000))]
    (and (= "INDEL" (:type vc))
         (splittable-size? vc)
         (or (has-multiple-nonref-alleles? vc)
             (has-ref-padding-mismatch? vc)))))

(defn- contains-indel? [alleles i]
  (when (< i (count (first alleles)))
    (contains? (set (map #(str (nth % i)) alleles)) "-")))

(defn- starts-an-indel? [alleles i]
  (contains-indel? alleles (inc i)))

(defn- gap-end? [alleles i]
  (and (pos? i)
       (not (contains-indel? alleles i))
       (contains-indel? alleles (dec i))))

(defn- is-match? [alleles i]
  (= 1 (count (set (map #(str (nth % i)) alleles)))))

(defn- split-alleles
  "Detect single call SNP variants within a MNP genotype.
  Handles reference no-variant padding bases on the 5' end of
  the sequence, writing only variants at the adjusted positions."
  [vc alleles & {:keys [prev-pad]}]
  (letfn [(is-internal-indel? [alleles i]
            (and (pos? i)
                 (contains-indel? alleles i)
                 (is-match? alleles (dec i))))
          (is-anchor-mismatch? [alleles i]
            (and (= 1 i)
                 (not (is-match? alleles i))
                 (is-match? alleles 0)))
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
          (has-nopad-five-indel? [alleles i]
            (and (is-fiveprime-indel? alleles i)
                 (contains-indel? alleles i)))
          (extend-indels [alleles i]
            {:start (if (or (is-internal-indel? alleles i)
                            (is-fiveprime-indel? alleles i))
                      (max (dec i) 0)
                      i)
             :end (inc (or (last (take-while #(or (contains-indel? alleles %)
                                                  (starts-an-indel? alleles %))
                                             (range i (count (first alleles)))))
                           i))})
          (ref-and-alt-alleles [cur-alleles]
            (let [refa (first cur-alleles)
                  alts (map (fn [x]
                              (if (= (.getDisplayString x) (.getDisplayString refa))
                                refa x))
                            (rest cur-alleles))]
              {:ref refa :alts alts}))
          (extract-variants [alleles pos]
            (let [{:keys [start end]} (extend-indels alleles pos)
                  str-alleles (map #(-> (str (if (has-nopad-five-indel? alleles start) prev-pad "")
                                             (subs % start end))
                                        (string/replace "-" ""))
                                   alleles)
                  cur-alleles (map-indexed (fn [i x] (Allele/create x (= 0 i)))
                                           str-alleles)
                  size (let [base (.length (first cur-alleles))]
                         (if (some empty? str-alleles) base (dec base)))
                  w-gap-start (-> (first alleles) (subs 0 start) (string/replace "-" "") count)
                  ready-alleles (ref-and-alt-alleles cur-alleles)]
              {:offset (+ w-gap-start (if (has-nopad-five-indel? alleles start) -1 0))
               :end (+ w-gap-start size)
               :next-start end
               :size size
               :orig-alleles alleles
               :ref-allele (:ref ready-alleles)
               :alleles (:alts ready-alleles)}))]
    (remove nil?
            (loop [i 0 final []]
              (cond
               (>= i (-> alleles first count)) final
               (has-variant-base? alleles i)
               (let [next-var (extract-variants alleles i)]
                 (recur (:next-start next-var) (conj final next-var)))
               :else (recur (inc i) final))))))

(defn- genotype-w-alleles
  "Retrieve a new set of genotypes with the given alleles.
   Update genotypes from the VariantContext, copying the existing
   genotype and substituting in the provided alleles and phasing information."
  [vc alleles orig-alleles is-phased]
  (letfn [(get-new-allele [new-alleles orig-alleles]
            (let [old-map (into {} (map-indexed
                                    (fn [i x]
                                      [(string/replace x "-" "") i])
                                    orig-alleles))]
              (fn [old-allele]
                (if (.isNoCall old-allele)
                  old-allele
                  (nth new-alleles (get old-map (.getDisplayString old-allele)))))))
          (add-new-genotype [allele-mapper context genotype]
            (doto context
              (.replace (-> (GenotypeBuilder. genotype)
                            (.alleles (map allele-mapper (.getAlleles genotype)))
                            (.phased (or (.isPhased genotype) is-phased))
                            .make))))]
    (reduce (partial add-new-genotype (get-new-allele alleles orig-alleles))
            (-> vc .getGenotypes GenotypesContext/copy)
            (.getGenotypes vc))))

(defn- new-split-vc
  "Create a new VariantContext as a subset of an existing variant.
   `allele-info` specifies the location size and alleles for the new variant:
   `{:offset :size :ref-allele :alleles}`"
  [vc i allele-info]
  (let [pos (+ (:offset allele-info) (.getStart vc))
        all-alleles (cons (:ref-allele allele-info) (:alleles allele-info))]
    (-> (VariantContextBuilder. vc)
        (.start pos)
        (.stop (+ pos (get allele-info :size 0)))
        (.genotypes (genotype-w-alleles vc all-alleles (:orig-alleles allele-info)
                                        (> i 0)))
        (.alleles (set all-alleles))
        (.make))))

(defn- split-mnp
  "Split a MNP into individual alleles"
  [vc]
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
            (let [n (count seqs)
                  gap-free (remove all-gap? (partition n (apply interleave (take n seqs))))]
              (map (fn [i]
                     (string/join "" (map #(nth % i) gap-free)))
                   (range n))))]
    (let [align-args (to-array [(SimpleGapPenalty. 20 1)])
          base-align (map #(.getSequenceAsString %)
                          (-> (map #(DNASequence. %) seqs)
                              (Alignments/getMultipleSequenceAlignment align-args)
                              .getAlignedSequences))
          orig-align (sort-by (original-seq-position seqs) (unique-aligns base-align))]
      (finalize-alignment orig-align))))

(defn- fix-gap-start-mismatch
  "Left align variants that start with a gap mismatch."
  [alleles]
  (letfn [(make-5-gap-wref [x]
            (let [anchor (subs x 0 1)
                  nogap-x (string/replace (subs x 1) "-" "")]
              (string/join "" (conj (vec (cons anchor
                                               (repeat (dec (- (count x) (count nogap-x))) "-")))
                                    nogap-x))))]
    (if (.contains (second alleles) "-")
      [(first alleles) (make-5-gap-wref (second alleles))]
      [(make-5-gap-wref (first alleles)) (second alleles)])))

(defn- left-align-complex
  "Ensure reference alignment gaps next to variants are consistently left aligned.
   Adjacent SNP and indels can have the SNP placed anywhere within the indel. This left
   aligns them to maintain anchoring via the 5' reference.
    ATCT  => ATCT
    AC--     A--C"
  [alleles]
  {:pre [(= 2 (count alleles))]
   :post [(= (count (first alleles))
             (count (first %)))]}
  (letfn [(gap-start-mismatch? [alleles i]
            (or (and (starts-an-indel? alleles i)
                     (not (is-match? alleles i))
                     (not (contains-indel? alleles i)))
                false))
          (gap-allele-type [alleles i]
            (cond
             (gap-start-mismatch? alleles i) :gs-mismatch
             (contains-indel? alleles i) :gap
             (gap-end? alleles i) :gap-end
             (is-match? alleles i) :match
             :else :mismatch))
          (split-at-match-gaps [[ann _]]
            (if (= :gap-end ann) ann :match))
          (get-region-allele [xs allele]
            (apply str (map #(nth allele (second %)) xs)))
          (get-region-alleles [alleles xs]
            (let [orig-alleles (map (partial get-region-allele xs) alleles)]
              (if (contains? (set (map first xs)) :gs-mismatch) 
                (fix-gap-start-mismatch orig-alleles)
                orig-alleles)))
          (concat-results [allele-parts]
            (vec (map #(apply str (map % allele-parts)) [first second])))]
    (concat-results
     (->> (map (fn [x] [(gap-allele-type alleles x) x]) (range (count (first alleles))))
          (partition-by split-at-match-gaps)
          (map (partial get-region-alleles alleles))))))

(defn- sanity-check-split-vcs
  "Confirm that new variants match correctly back to original.
   Catch any potential errors in splitting by ensuring reference coordinates
   and sequences match original."
  [vc new-vcs]
  (letfn [(get-vc-info [vc]
            (let [alleles (map #(.getDisplayString %) (.getAlleles vc))]
              {:start (.getStart vc)
               :alleles alleles}))
          (get-check-ref [orig new]
            (let [int-pos (- (:start new) (:start orig))
                  check-ref (first (:alleles new))]
              (if (neg? int-pos)
                [0 (subs check-ref (Math/abs int-pos))]
                [int-pos check-ref])))
          (check-split-vc [orig new]
            (let [[int-pos check-ref] (get-check-ref orig new)]
              (when (or (>= int-pos (count (first (:alleles orig))))
                        (neg? int-pos)
                        (not= (subs (first (:alleles orig)) int-pos (+ int-pos (count check-ref)))
                           check-ref))
                (throw (Exception. (format "Matching problem with split alleles: %s %s %s %s"
                                           (:chr vc) (:start vc) orig new))))))]
    (doall (map (partial check-split-vc (get-vc-info (:vc vc)))
                (map get-vc-info new-vcs)))))

(defn- split-complex-indel
  "Split complex indels into individual variant components."
  [vc ref]
  (let [prev-pad (or (extract-sequence ref (:chr vc) (dec (:start vc)) (dec (:start vc))) "N")
        ; Do not use reference sequence. Causes more trouble than aligning ref/alt directly.
        ref-seq nil ;(extract-sequence ref (:chr vc) (:start vc) (:end vc))
        alleles (split-alleles vc (->> (conj (get-vc-alleles vc) ref-seq)
                                       (remove empty?)
                                       (remove nil?)
                                       multiple-alignment
                                       (#(if (> (count %) 2) % (left-align-complex %))))
                               :prev-pad prev-pad)]
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
  (letfn [(strip-indel [vc i alleles]
            (let [start-pos (dec i)
                  ref-allele (subs (first alleles) start-pos)
                  cur-alleles (map #(Allele/create (subs % start-pos)
                                                   (= ref-allele (subs % start-pos)))
                                   alleles)]
              (new-split-vc vc 0 {:offset start-pos
                                  :size (dec (count ref-allele))
                                  :orig-alleles alleles
                                  :ref-allele (first cur-alleles)
                                  :alleles (rest cur-alleles)})))
          (variant-allele-pos [input-alleles]
            (let [str-alleles (map #(.getDisplayString %) input-alleles)
                  first-var-i (first (filter #(has-variant-base? str-alleles %)
                                             (range (apply max (map count str-alleles)))))]
              [str-alleles first-var-i]))
          (used-alt-alleles [vc]
            (let [genotype-alleles (set (mapcat :alleles (:genotypes vc)))]
              (filter #(contains? genotype-alleles %) (:alt-alleles vc))))]
    (let [alt-alleles (used-alt-alleles vc)
          [orig-alleles first-var-i] (variant-allele-pos (cons (:ref-allele vc)
                                                               alt-alleles))
          [_ nocall-i] (variant-allele-pos (cons (:ref-allele vc) alt-alleles))]
      (if (or (nil? first-var-i) (<= first-var-i 1)
              (nil? nocall-i) (<= nocall-i 1))
        (:vc vc)
        (strip-indel (:vc vc) first-var-i orig-alleles)))))

;; ## VCF file conversion
;; Process entire files, normalizing complex variations

(defn- round-vc-start
  "Round the start value of a VC to the nearest ten million.
   This heuristic will cause problems with out of order variant
   contexts that span this split junction (9999999 and 10000000) but
   saves having to work with overlapping groups allowing streaming."
  [vc]
  (let [rounder 10000000.0]
    {:chr (.getChr vc)
     :pos (-> (.getStart vc)
              (/ rounder)
              Math/floor
              (* rounder)
              int)}))

(defn- sort-vc-group
  "Sort a group of variant contexts by start position.
   Ensures that post-normalized variant contexts
   sort correctly within the blocks defined by round-vc-start."
  [vcs]
  (sort-by #(.getStart %) vcs))

(defn- get-normalized-vcs
  "Lazy list of variant context with MNPs split into single genotypes and indels stripped."
  [vc-iter ref]
  (letfn [(process-vc [vc]
            (condp = (:type vc)
              "MNP" (split-mnp vc)
              "INDEL" (if (is-multi-indel? vc)
                        (split-complex-indel vc ref)
                        [(maybe-strip-indel vc)])
              [(:vc vc)]))]
    (->> (mapcat process-vc vc-iter)
         (partition-by round-vc-start)
         (mapcat sort-vc-group))))

(defn left-align-variants
  "Left align variants in an input VCF file for a standard representation.
  Checks final line count of prepared file, returning left-aligned files
  only when converting every variant in the input."
  [in-file ref & {:keys [out-dir rerun?]}]
  (letfn [(line-count [f]
            (with-open [rdr (reader f)]
              (count (remove #(.startsWith % "#") (line-seq rdr)))))]
    (let [file-info {:out-vcf (fsp/add-file-part in-file "leftalign" out-dir)}
          args ["-R" ref "-o" :out-vcf "--variant" in-file]]
      (when (and rerun? (fs/exists? (:out-vcf file-info)))
        (fs/delete (:out-vcf file-info)))
      (broad/run-gatk "LeftAlignVariants" args file-info {:out [:out-vcf]})
      (if (= (line-count in-file) (line-count (:out-vcf file-info)))
        (:out-vcf file-info)
        in-file))))

(defn normalize-variants
  "Convert MNPs and indels into normalized representation."
  ([in-file ref]
     (normalize-variants in-file ref nil))
  ([in-file ref out-dir & {:keys [out-fname]}]
     (let [base-name (if (nil? out-fname) (fsp/remove-zip-ext in-file) out-fname)
           out-file (fsp/add-file-part base-name "nomnp" out-dir)
           out-pre-file (fsp/add-file-part base-name "worknomnp" out-dir)]
       (when (itx/needs-run? out-file)
         (when (fs/exists? out-pre-file)
           (fs/delete out-pre-file))
         (let [la-file (left-align-variants in-file ref :out-dir out-dir :rerun? true)]
           (with-open [vcf-iter (get-vcf-iterator la-file ref)]
             (write-vcf-w-template in-file {:out out-pre-file}
                                   (get-normalized-vcs (parse-vcf vcf-iter) ref)
                                   ref))
           (fs/rename (left-align-variants out-pre-file ref :out-dir out-dir :rerun? true)
                      out-file)))
       out-file)))
