(ns bcbio.variation.annotate.entropy
  "Calculate Shannon entropy for flanking sequence surrounding variants.
   Used to identify low-complexity repeat regions in variants.
   Based on 'vcfentropy' from Erik Garrison's vcflib:
   https://github.com/ekg/vcflib"
  (:import [org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation]
           [org.broadinstitute.variant.vcf VCFInfoHeaderLine VCFHeaderLineType])
  (:gen-class
   :name bcbio.variation.annotate.entropy.ShannonEntropy
   :extends org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation))


;; ## Shannon entropy
;; From John Lawrence Aspden's information theory posts
;; https://github.com/johnlawrenceaspden/hobby-code/blob/master/averygentleintroduction-part5.clj

(defn bits [n]
  "How many bits to represent n alternatives? Fractions allowed! Also know as log2."
  (/ (Math/log n) (Math/log 2)))

(defn shannon-entropy [P]
  (let [odds (map second P)
        total (reduce + odds)
        bits-aliased (/ (reduce + (map * odds (map bits odds))) total)]
    (- (bits total) bits-aliased)))

(defn seq-entropy
  "Calculate entropy of a sequence based on distribution of dimers.
   Splits sequence into all dimer 2bp windows, calculates frequency
   of each dimer and then feeds distribution to shannon calculation."
  [seq]
  (->> seq
       (partition 2 1)
       frequencies
       shannon-entropy))

;; ## Helper function

(defn get-flank-seq
  "Retrieve sequence surrounding the current variant, with nbp flanking sequence."
  [ref-context nbp]
  (letfn [(subset-region [x]
            (let [want-size (inc (* 2 nbp))
                  end-subtract (/ (- (count x) want-size) 2)]
              (subs x end-subtract (- (count x) end-subtract))))]
    (->> ref-context
         .getBases
         (map char)
         (apply str)
         subset-region)))

;; ## GATK walker

(def flank-bp 12)

(defn -getKeyNames [_]
  ["Entropy"])

(defn -getDescriptions [_]
  [(VCFInfoHeaderLine. "Entropy" 1 VCFHeaderLineType/Float
                       (format "Shannon entropy of variant flanking regions, %sbp on both sides"
                               flank-bp))])

(defn -annotate
  "Retrieve flanking region surrounding variant and calculate entropy."
  [_ _ _ ref _ _ _]
  {"Entropy" (->> (get-flank-seq ref flank-bp)
                  seq-entropy
                  (format "%.2f"))})
