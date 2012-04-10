(ns bcbio.variation.structural
  "Handle structural variations for larger insertions, deletions and
  genome rearrangements."
  (:import [org.broadinstitute.sting.utils.codecs.vcf VCFCodec]
           [org.broadinstitute.sting.utils.variantcontext VariantContextBuilder
            Allele]
           [java.io StringReader]
           [net.sf.picard.util IntervalTree])
  (:use [bcbio.variation.variantcontext :only [get-vcf-source parse-vcf
                                               from-vc write-vcf-w-template]]
        [bcbio.variation.callable :only [get-bed-source]])
  (:require [clojure.string :as string]
            [fs.core :as fs]
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

(defn value-from-attr
  "Retrieve normalized integer values from an attribute."
  [vc attr-name attr-index]
  (-> vc
      :attributes
      (get attr-name (repeat (inc attr-index) "0"))
      (#(if (string? %) [%] %))
      (nth attr-index)
      (Integer/parseInt)
      Math/abs))

(defn parse-vcf-sv
  "Parse VCF file returning structural variants with confidence intervals.
  The :out-format keyword specifies how to return the parsed structural variants:
   - :itree -- Interval tree for variant lookup by chromosome and start/end.
   - default -- List of variants (non-lazy)."
  [vcf-file ref-file & {:keys [out-format interval-file]}]
  (letfn [(start-adjust [vc]
            (value-from-attr vc "CIPOS" 0))
          (end-adjust [vc]
            (value-from-attr vc "CIEND" 1))
          (updated-sv-vc [cur-vc]
            (when-let [sv-type (get-sv-type cur-vc)]
              (-> cur-vc
                  (assoc :start-ci (- (:start cur-vc) (start-adjust cur-vc)))
                  (assoc :end-ci (+ (:end cur-vc) (end-adjust cur-vc)))
                  (assoc :sv-type sv-type))))
          (prep-itree [vc-iter]
            (reduce (fn [coll vc]
                      (assoc coll (:chr vc)
                             (doto (get coll (:chr vc) (IntervalTree.))
                               (.put (:start-ci vc) (inc (:end-ci vc)) vc))))
                    {} vc-iter))
          (in-intervals? [bed-source vc]
            (or (instance? StringReader bed-source)
                (not (nil? (first (.query bed-source (:chr vc) (:start-ci vc) (:end-ci vc)))))))]
    (with-open [vcf-source (get-vcf-source vcf-file ref-file :codec (structural-vcfcodec))
                bed-source (if-not (nil? interval-file) (get-bed-source interval-file)
                                   (StringReader. ""))]
      (let [vs-iter (filter (partial in-intervals? bed-source)
                            (keep updated-sv-vc (parse-vcf vcf-source)))]
        (case out-format
          :itree (prep-itree vs-iter)
          (vec vs-iter))))))

;; ## Concordance checking

(defn- coords-with-ci
  "Retrieve the coordinates of a variant context with confidence intervals."
  [s e vc]
  (letfn [(get-ci-range [orig attr]
            [(- orig (value-from-attr vc attr 0))
             (+ orig (value-from-attr vc attr 1))])]
    [(get-ci-range s "CIPOS")
     (get-ci-range e "CIEND")]))

(defn sv-ends-overlap?
  "Check if coordinates from two structural variants overlap.
  Considered an overlap if the two confidence intervals
  have shared bases."
  [[[s1 e1] [s2 e2]]]
  (or (and (>= s2 s1) (<= s2 e1))
      (and (>= e2 s1) (<= e2 e1))))

(defn sv-len-concordant?
  "Check for concordance of variants based on reported length: deletions and inversions."
  [sv1 sv2]
  (letfn [(get-start-end [x]
            (let [start (:start x)
                  end (max (+ start (value-from-attr x "SVLEN" 0))
                           (:end x))]
              (coords-with-ci start end x)))]
    (every? sv-ends-overlap?
            (partition 2 (interleave (get-start-end sv1) (get-start-end sv2))))))

(defn sv-concordant?
  "Check if structural variants are concordant."
  [sv1 sv2]
  (and (apply = (map :sv-type [sv1 sv2]))
       (case (:sv-type sv1)
         :DEL (sv-len-concordant? sv1 sv2)
         :INS true
         :INV (sv-len-concordant? sv1 sv2)
         :DUP true
         :BND true
         (throw (Exception. (str "Structural variant type not handled: " (:sv-type sv1)))))))

(defn get-itree-overlap
  "Lazy sequence of items that overlap a region in a nested IntervalTree."
  [itree chrom start end]
  (letfn [(itree-seq [iter]
            (lazy-seq
             (when (.hasNext iter)
               (cons (.getValue (.next iter)) (itree-seq iter)))))]
    (-> itree
        (get chrom)
        (.overlappers start end)
        itree-seq)))

(defn- find-concordant-svs
  "Compare two structural variant files, returning variant contexts keyed by concordance."
  [fname1 fname2 ref interval-file]
  (let [cmp-tree (parse-vcf-sv fname2 ref :out-format :itree :interval-file interval-file)]
    (letfn [(check-sv-concordance [vc]
              (let [matches (filter (partial sv-concordant? vc)
                                    (get-itree-overlap cmp-tree (:chr vc)
                                                       (:start-ci vc) (inc (:end-ci vc))))]
              [(if (seq matches) :concordant :discordant1) (:vc vc)]))]
      (map check-sv-concordance (parse-vcf-sv fname1 ref :interval-file interval-file)))))

(defn compare-sv
  "Compare structural variants, producing concordant and discordant outputs"
  [sample c1 c2 ref & {:keys [out-dir interval-file]}]
  (let [base-out (str (fs/file (if (nil? out-dir) (fs/parent (:file c1)) out-dir)
                               (str sample "-%s-%s-%s.vcf")))
        out-files {:concordant (format base-out (:name c1) (:name c2) "svconcordance")
                   :discordant1 (format base-out (:name c1) (:name c2) "svdiscordance")
                   :discordant2 (format base-out (:name c2) (:name c1) "svdiscordance")}]
    (when (itx/needs-run? (vals out-files))
      (write-vcf-w-template (:file c1) out-files
                            (find-concordant-svs (:file c1) (:file c2) ref interval-file)
                            ref))
    out-files))
