(ns bcbio.variation.structural
  "Handle structural variations for larger insertions, deletions and
  genome rearrangements."
  (:import [org.broadinstitute.sting.utils.codecs.vcf VCFCodec]
           [org.broadinstitute.sting.utils.variantcontext VariantContextBuilder
            Allele])
  (:use [bcbio.variation.variantcontext :only [get-vcf-source parse-vcf
                                               from-vc]])
  (:require [clojure.string :as string]
            [bcbio.run.itx :as itx]))

(defn get-sv-type
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

(defn nochange-alt?
  "Check a VCF input line for identical REF and ALT calls"
  [line]
  (let [parts (string/split line #"\t")]
    (= (nth parts 3) (nth parts 4))))

(defn- check-sv-line
  "Check SV inputs for validity, fixing or filtering where possible.
  Fixes:
    - identical ref/alt calls: an apparent SV no-call"
  [line]
  (cond
   (.startsWith line "#") line
   (nochange-alt? line) nil
   :else line))

(defn structural-vcfcodec []
  "Provide VCFCodec decoder that returns structural variants expanded
  to include confidence regions."
  (letfn [(pos-from-attr [vc attr-name attr-index]
            (-> vc
                :attributes
                (get attr-name ["0" "0"])
                (nth attr-index)
                (Integer/parseInt)
                Math/abs))
          (start-adjust [vc]
            (pos-from-attr vc "CIPOS" 0))
          (end-adjust [vc]
            (pos-from-attr vc "CIEND" 1))
          (fix-ref [allele new-bases]
            (Allele/create (apply str (cons (.getBaseString allele) (repeat new-bases "N")))
                           true))
          (update-pos [vc sv-type]
            (let [start-pad (start-adjust vc)
                  end-pad (end-adjust vc)]
              (-> (VariantContextBuilder. (:vc vc))
                  (.start (- (:start vc) start-pad))
                  (.stop (+ (:end vc) end-pad))
                  (.alleles (set (cons (fix-ref (:ref-allele vc) (+ start-pad end-pad))
                                       (:alt-alleles vc))))
                  .make)))
          (updated-sv-vc [vc]
            (let [cur-vc (from-vc vc)]
              (when-let [sv-type (get-sv-type cur-vc)]
                (update-pos cur-vc sv-type))))]
    (proxy [VCFCodec] []
      (decode [line]
        (when-let [work-line (check-sv-line line)]
          (when-let [vc (proxy-super decode work-line)]
            (updated-sv-vc vc))))
      (decodeLoc [line]
        (when-let [work-line (check-sv-line line)]
          (when-let [vc (proxy-super decode work-line)]
            (updated-sv-vc vc)))))))

(defn parse-sv-vcf [vcf-file ref-file]
  (itx/remove-path (str vcf-file ".idx"))
  (with-open [vcf-source (get-vcf-source vcf-file ref-file :codec (structural-vcfcodec))]
    (doseq [vc (parse-vcf vcf-source)]
      (println vc))))
