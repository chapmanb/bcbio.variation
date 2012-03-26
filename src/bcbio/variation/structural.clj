(ns bcbio.variation.structural
  "Handle structural variations for larger insertions, deletions and
  genome rearrangements."
  (:import [org.broadinstitute.sting.utils.codecs.vcf VCFCodec]
           [org.broadinstitute.sting.utils.variantcontext VariantContextBuilder
            Allele])
  (:use [bcbio.variation.variantcontext :only [get-vcf-source parse-vcf
                                               from-vc]])
  (:require [clojure.string :as string]
            [bcbio.run.itx :as itx]
            [bcbio.itree :as itree]))

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

(defn structural-vcfcodec []
  "Provide VCFCodec decoder that returns structural variants.
  Check SV inputs for validity, fixing or filtering where possible.
  Fixes:
    - identical ref/alt calls: an apparent SV no-call"
  (letfn [(check-sv-line [line]
            (cond
             (.startsWith line "#") line
             (nochange-alt? line) nil
             :else line))]
    (proxy [VCFCodec] []
      (decode [line]
        (when-let [work-line (check-sv-line line)]
          (when-let [vc (proxy-super decode work-line)]
            vc)))
      (decodeLoc [line]
        (when-let [work-line (check-sv-line line)]
          (when-let [vc (proxy-super decode work-line)]
            vc))))))

(defn parse-vcf-sv
  "Parse VCF file returning structural variants with confidence intervals.
  The :out-format keyword specifies how to return the parsed structural variants:
   - :itree -- Interval tree for variant lookup by chromosome and start/end.
   - default -- List of variants (non-lazy)."
  [vcf-file ref-file & {:keys [out-format]}]
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
          (updated-sv-vc [cur-vc]
            (when-let [sv-type (get-sv-type cur-vc)]
              (-> cur-vc
                  (assoc :start-ci (- (:start cur-vc) (start-adjust cur-vc)))
                  (assoc :end-ci (+ (:end cur-vc) (end-adjust cur-vc))))))
          (prep-itree [vc-iter]
            (reduce (fn [coll vc] (assoc coll (:chr vc)
                                         (itree/iassoc (get coll (:chr vc) itree/empty-interval-map)
                                                       (:start-ci vc) (:end-ci vc) vc)))
                    {} vc-iter))]
    (with-open [vcf-source (get-vcf-source vcf-file ref-file :codec (structural-vcfcodec))]
      (let [vs-iter (keep updated-sv-vc (parse-vcf vcf-source))]
        (case out-format
          :itree (prep-itree vs-iter)
          (vec vs-iter))))))
