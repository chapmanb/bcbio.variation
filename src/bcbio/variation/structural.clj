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
        [bcbio.variation.callable :only [get-bed-source]]
        [ordered.map :only [ordered-map]])
  (:require [clojure.string :as string]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]))

;; ## Interval tree lookup

(defn prep-itree
  "Retrieve an Interval with the specified start/end keywords."
  [vc-iter start-kw end-kw]
  (reduce (fn [coll vc]
            (assoc coll (:chr vc)
                   (doto (get coll (:chr vc) (IntervalTree.))
                     (.put (get vc start-kw) (inc (get vc end-kw)) vc))))
          (ordered-map) vc-iter))

(defn- itree-seq
  "Convert IntervalTree Iterator into clojure seq."
  [iter]
  (lazy-seq
   (when (.hasNext iter)
     (cons (.getValue (.next iter)) (itree-seq iter)))))

(defn get-itree-overlap
  "Lazy sequence of items that overlap a region in a nested IntervalTree."
  [itree chrom start end]
  (-> itree
      (get chrom)
      (.overlappers start end)
      itree-seq))

(defn get-itree-all
  "Lazy sequence of all items in an IntervalTree."
  [itree]
  (flatten
   (for [item (vals itree)]
     (itree-seq (.iterator item)))))

(defn remove-itree-vc
  "Remove variant context from an IntervalTree"
  [itree chr start end]
  (if (not-any? nil? [chr start end])
    (assoc itree chr
           (doto (get itree chr)
             (.remove start (inc end))))
    itree))

;; ## Structural variation helpers

(defn get-sv-type
  "Determine the type of a structural variant. Expected types are:

    - DEL: Deletion
    - INS: Insertion
    - DUP: Duplication
    - INV: Inversion
    - BND: Breakpoint end; paired with second variant
    - CNV: Copy number variation
    - nil: Not a structural variant."
  [vc params]
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
          (> (max-allele-size vc) (get params :min-indel 100))) (indel-type vc)
     (= "SYMBOLIC" (:type vc)) (alt-sv-type vc)
     :else nil)))

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

;; ## Parsing structural variants

(defn value-from-attr
  "Retrieve normalized integer values from an attribute."
  ([vc attr-name]
     (value-from-attr vc attr-name 0))
  ([vc attr-name attr-index]
      (-> vc
          :attributes
          (get attr-name (repeat (inc attr-index) "0"))
          (#(if (string? %) [%] %))
          (nth attr-index)
          (Integer/parseInt)
          Math/abs)))

(defn parse-vcf-sv
  "Parse VCF file returning structural variants with confidence intervals.
  The :out-format keyword specifies how to return the parsed structural variants:
   - :itree -- Interval tree for variant lookup by chromosome and start/end.
   - default -- List of variants (non-lazy)."
  [vcf-file ref-file & {:keys [out-format interval-file params]
                        :or {params {}}}]
  (letfn [(start-adjust [vc]
            (value-from-attr vc "CIPOS" 0))
          (end-adjust [vc]
            (value-from-attr vc "CIEND" 1))
          (updated-sv-vc [cur-vc]
            (when-let [sv-type (get-sv-type cur-vc params)]
              (-> cur-vc
                  (assoc :start-ci (- (:start cur-vc) (start-adjust cur-vc)))
                  (assoc :end-ci (+ (:end cur-vc) (end-adjust cur-vc)))
                  (assoc :sv-type sv-type))))
          (in-intervals? [bed-source vc]
            (or (instance? StringReader bed-source)
                (not (nil? (first (.query bed-source (:chr vc) (:start-ci vc) (:end-ci vc)))))))]
    (with-open [vcf-source (get-vcf-source vcf-file ref-file :codec (structural-vcfcodec))
                bed-source (if-not (nil? interval-file) (get-bed-source interval-file)
                                   (StringReader. ""))]
      (let [vs-iter (filter (partial in-intervals? bed-source)
                            (keep updated-sv-vc (parse-vcf vcf-source)))]
        (case out-format
          :itree (prep-itree vs-iter :start-ci :end-ci)
          (vec vs-iter))))))

;; ## Concordance checking

(defn get-ci-start-end
  "Retrieve start and end with confidence intervals for a variation.
  length-fn returns the length of the item or a string for well-known items."
  [vc length-fn]
  (letfn [(get-ci-range [orig attr]
            [(- orig (value-from-attr vc attr 0))
             (+ orig (value-from-attr vc attr 1))])]
    (let [start (:start vc)
          length (length-fn vc)
          end (if (string? length) length
                  (max (+ start length) (:end vc)))]
      [(get-ci-range start "CIPOS")
       (if (string? end) end
           (get-ci-range end "CIEND"))])))

(defmulti sv-ends-overlap?
  "Check if coordinates from two structural variants overlap.
  Considered an overlap if the two confidence intervals
  have shared bases."
  (fn [[end1 end2]] (type end1)))

(defmethod sv-ends-overlap? clojure.lang.PersistentVector
  [[[s1 e1] [s2 e2]]]
  (or (and (>= s2 s1) (<= s2 e1))
      (and (>= e2 s1) (<= e2 e1))))

(defmethod sv-ends-overlap? java.lang.String
  [[end1 end2]]
  (= end1 end2))

(defn- length-from-svlen [x] (value-from-attr x "SVLEN"))

(defn- insertion-length
  "Length of insertion variation, handling ALT allele, INSEQ
  and well-known named insertions."
  [x]
  (letfn [(get-insseq [x]
            (-> x :attributes (get "INSSEQ")))
          (length-by-insert-name [alt-allele]
            (cond
             (.startsWith alt-allele "<INS:ME:") (-> alt-allele
                                                     (subs 1 (dec (count alt-allele)))
                                                     (string/split #":")
                                                     last)
             :else (throw (Exception. (str "Unknown insert allele" alt-allele)))))
          (get-allele-insert [x]
            (let [alt-allele (-> x :alt-alleles first .getDisplayString)]
              (if (.startsWith alt-allele "<")
                (length-by-insert-name alt-allele)
                (dec (count alt-allele)))))]
    (if-let [seq (get-insseq x)]
      (count seq)
      (get-allele-insert x))))

(defn- duplication-length
  "Length of duplication variation, handling SVLEN and END."
  [x]
  (max (length-from-svlen x)
       (- (value-from-attr x "END") (:start x))))

(defn sv-len-concordant?
  "Check for concordance of variants based on reported length:
  handles deletions, inversions. insertions and duplications.
  length-fn is a custom function to retrieve the variation length."
  [sv1 sv2 length-fn]
  (letfn []
    (every? sv-ends-overlap?
            (partition 2 (interleave (get-ci-start-end sv1 length-fn)
                                     (get-ci-start-end sv2 length-fn))))))

(defn sv-partial-match?
  "Allow partial matching of structural variants based on overlaps"
  [sv1 sv2]
  (and (= (:chr sv1) (:chr sv2)))
  (println (map (juxt :start :end) [sv1 sv2])))

(defn sv-concordant?
  "Check if structural variants are concordant."
  [params sv1 sv2]
  (and (apply = (map :sv-type [sv1 sv2]))
       (case (keyword (get params :method "end-ci"))
         :end-ci (case (:sv-type sv1)
                   :DEL (sv-len-concordant? sv1 sv2 length-from-svlen)
                   :INS (sv-len-concordant? sv1 sv2 insertion-length)
                   :INV (sv-len-concordant? sv1 sv2 length-from-svlen)
                   :DUP (sv-len-concordant? sv1 sv2 duplication-length)
                   :BND false
                   (throw (Exception. (str "Structural variant type not handled: "
                                           (:sv-type sv1)))))
         :partial-overlap (case (:sv-type sv1)
                    (:DEL :INS) (sv-partial-match? sv1 sv2)
                    (throw (Exception. (str "Structural variant type not handled: "
                                           (:sv-type sv1))))))))

(defn- find-concordant-svs
  "Compare two structural variant files, returning variant contexts keyed by concordance."
  [fname1 fname2 ref interval-file params]
  (let [cmp-tree (atom (parse-vcf-sv fname2 ref :out-format :itree :interval-file interval-file
                                     :params params))]
    (letfn [(check-sv-concordance [vc]
              (let [matches (filter (partial sv-concordant? params vc)
                                    (get-itree-overlap @cmp-tree (:chr vc)
                                                       (:start-ci vc) (inc (:end-ci vc))))]
                (doseq [m-vc matches]
                  (reset! cmp-tree (remove-itree-vc @cmp-tree (:chr m-vc)
                                                    (:start m-vc) (:end m-vc))))
                [(if (seq matches) :concordant :discordant1) (:vc vc)]))
            (remaining-cmp-svs [itree]
              (partition 2
                         (interleave (repeat :discordant2) (map :vc (get-itree-all itree)))))]

      (concat
       (map check-sv-concordance (parse-vcf-sv fname1 ref :interval-file interval-file
                                               :params params))
       (remaining-cmp-svs @cmp-tree)))))

(defn compare-sv
  "Compare structural variants, producing concordant and discordant outputs"
  [sample c1 c2 ref & {:keys [out-dir interval-file params]
                       :or {params {}}}]
  (let [base-out (str (fs/file (if (nil? out-dir) (fs/parent (:file c1)) out-dir)
                               (str sample "-%s-%s-%s.vcf")))
        out-files {:concordant (format base-out (:name c1) (:name c2) "svconcordance")
                   :discordant1 (format base-out (:name c1) (:name c2) "svdiscordance")
                   :discordant2 (format base-out (:name c2) (:name c1) "svdiscordance")}]
    (when (itx/needs-run? (vals out-files))
      (write-vcf-w-template (:file c1) out-files
                            (find-concordant-svs (:file c1) (:file c2) ref interval-file
                                                 params)
                            ref))
    out-files))

;; ## Utility functions

(defn- find-overlapping-svs
  "Lazy stream of structural variants overlapping in both inputs."
  [f1 f2 ref params]
  (letfn [(find-overlaps [cmp-tree vc]
            (when-let [cmp-vc (first (get-itree-overlap cmp-tree (:chr vc)
                                                        (:start-ci vc) (inc (:end-ci vc))))]
              [:out1 (:vc vc) :out2 (:vc cmp-vc)]))]
    (let [cmp-tree (parse-vcf-sv f2 ref :out-format :itree :params params)]
      (->> (parse-vcf-sv f1 ref :params params)
           (map (partial find-overlaps cmp-tree))
           (remove nil?)
           flatten
           (partition 2)))))

(defn overlapping-svs
  "Prepare VCF files of only overlapping structural variants present in both."
  [f1 f2 ref params]
  (let [out-files {:out1 (itx/add-file-part f1 "overlap")
                   :out2 (itx/add-file-part f2 "overlap")}]
    (when (itx/needs-run? (vals out-files))
      (write-vcf-w-template f1 out-files (find-overlapping-svs f1 f2 ref params) ref))
    out-files))
