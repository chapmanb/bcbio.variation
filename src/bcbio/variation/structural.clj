(ns bcbio.variation.structural
  "Handle structural variations for larger insertions, deletions and
  genome rearrangements."
  (:import [org.broadinstitute.sting.utils.codecs.vcf VCFCodec]
           [org.broadinstitute.sting.utils.variantcontext VariantContextBuilder
            Allele]
           [java.io StringReader]
           [net.sf.picard.util IntervalTree])
  (:use [clojure.set :only [intersection]]
        [ordered.map :only [ordered-map]]
        [bcbio.variation.variantcontext :only [get-vcf-source parse-vcf
                                               from-vc write-vcf-w-template]]
        [bcbio.variation.callable :only [get-bed-source]]
        )
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
  "Convert IntervalTree Iterator into clojure seq.
  Catch deleted sequences and continue ignoring the deleted node."
  [iter]
  (lazy-seq
   (when (.hasNext iter)
     (try
       (cons (.getValue (.next iter)) (itree-seq iter))
       (catch java.util.ConcurrentModificationException e
         (itree-seq iter))))))

(defn get-itree-overlap
  "Lazy sequence of items that overlap a region in a nested IntervalTree."
  [itree chrom start end]
  (let [chr-itree (get itree chrom)]
    (if (nil? chr-itree)
      []
      (itree-seq (.overlappers chr-itree start end)))))

(defn get-itree-all
  "Lazy sequence of all items in an IntervalTree."
  [itree]
  (flatten
   (for [item (vals itree)]
     (sort-by :start
              (itree-seq (.iterator item))))))

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
            (if (> (.length (:ref-allele vc))
                   (apply max (map #(.length %) (:alt-alleles vc)))) :DEL :INS))
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
          (> (max-allele-size vc) (get params :min-indel 10))) (indel-type vc)
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

;; ## Concordance checking

(defmulti sv-ends-overlap?
  "Check if coordinates from two structural variants overlap.
  Considered an overlap if the two confidence intervals
  have shared bases."
  (fn [[end1 end2]] (type end1)))

(defmethod sv-ends-overlap? clojure.lang.PersistentVector
  [[[s1 e1] [s2 e2]]]
  (seq (intersection (set (range s1 (inc e1)))
                     (set (range s2 (inc e2))))))

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

(defn- deletion-length
  "Length of deletion variations, handling SVLEN and allele specifications."
  [vc]
  (let [svlen (length-from-svlen vc)]
    (if (pos? svlen)
      svlen
      (- (-> vc :ref-allele .length)
         (apply min (map #(.length %) (:alt-alleles vc)))))))

(defn- duplication-length
  "Length of duplication variation, handling SVLEN and END."
  [x]
  (max (length-from-svlen x)
       (- (value-from-attr x "END") (:start x))))

(defn- get-sv-length
  "Retrieve length of a structural variant for different variation types."
  [vc]
  (case (:sv-type vc)
         :DEL (deletion-length vc)
         :INS (insertion-length vc)
         :INV (length-from-svlen vc)
         :DUP (duplication-length vc)
         :CNV (duplication-length vc)
         :BND 0
         (throw (Exception. (str "Structural variant type not handled: "
                                 (:sv-type vc))))))

(defn- get-ci-start-end
  "Retrieve start and end with confidence intervals for a variation."
  [vc params & {:keys [allow-named?]}]
  (letfn [(get-ci-range [orig attr default-ci]
            (let [left-ci (value-from-attr vc attr 0)
                  right-ci (value-from-attr vc attr 1)]
              [(- orig (if (pos? left-ci) left-ci default-ci))
               (+ orig (if (pos? right-ci) right-ci default-ci))]))
          (get-default-ci [length]
            (let [default (if-let [x (-> (:default-cis params) first second)] x 0)
                  by-length (when-not (string? length)
                              (second (first (drop-while #(< (first %) length)
                                                         (:default-cis params)))))]
              (if-not (nil? by-length) by-length default)))]
    (let [start (:start vc)
          length (get-sv-length vc)
          default-ci (get-default-ci length)
          end (cond
               (and allow-named? (string? length)) length
               (string? length) (:end vc)
               :else (max (+ start length) (:end vc)))]
      [(get-ci-range start "CIPOS" default-ci)
       (if (string? end) end
           (get-ci-range end "CIEND" default-ci))])))

(defn- sv-len-concordant?
  "Check for concordance of variants based on reported length:
  handles deletions, inversions. insertions and duplications."
  [sv1 sv2 params]
  (every? sv-ends-overlap?
          (partition 2 (interleave (get-ci-start-end sv1 params :allow-named? true)
                                   (get-ci-start-end sv2 params :allow-named? true)))))

(defn sv-concordant?
  "Check if structural variants are concordant."
  [params sv1 sv2]
  (and (apply = (map :sv-type [sv1 sv2]))
       (case (:sv-type sv1)
         (:DEL :INS :INV :DUP) (sv-len-concordant? sv1 sv2 params)
         :BND false
         (throw (Exception. (str "Structural variant type not handled: "
                                 (:sv-type sv1)))))))

;; ## Parsing structural variants

(defn parse-vcf-sv
  "Parse VCF file returning structural variants with confidence intervals.
  The :out-format keyword specifies how to return the parsed structural variants:
   - :itree -- Interval tree for variant lookup by chromosome and start/end.
   - default -- List of variants (non-lazy)."
  [vcf-file ref-file & {:keys [out-format interval-file params]
                        :or {params {}}}]
  (letfn [(updated-sv-vc [cur-vc]
            (when-let [sv-type (get-sv-type cur-vc params)]
              (let [[start-cis end-cis] (get-ci-start-end (assoc cur-vc :sv-type sv-type)
                                                          params)]
                (-> cur-vc
                    (assoc :start-ci (first start-cis))
                    (assoc :end-ci (second end-cis))
                    (assoc :sv-type sv-type)))))
          (in-intervals? [bed-source vc]
            (or (instance? StringReader bed-source)
                (not (nil? (first (.query bed-source (:chr vc) (:start-ci vc) (:end-ci vc)))))))]
    (with-open [vcf-source (get-vcf-source vcf-file ref-file :codec (structural-vcfcodec))
                bed-source (if-not (nil? interval-file) (get-bed-source interval-file ref-file)
                                   (StringReader. ""))]
      (let [vs-iter (filter (partial in-intervals? bed-source)
                            (keep updated-sv-vc (parse-vcf vcf-source)))]
        (case out-format
          :itree (prep-itree vs-iter :start-ci :end-ci)
          (vec vs-iter))))))

(defn- find-concordant-svs
  "Compare two structural variant files, returning variant contexts keyed by concordance."
  [fname1 fname2 disc-kwds ref interval-file params]
  (let [cmp-tree (atom (parse-vcf-sv fname2 ref :out-format :itree :interval-file interval-file
                                     :params params))]
    (letfn [(check-sv-concordance [vc]
              (let [matches (filter (partial sv-concordant? params vc)
                                    (get-itree-overlap @cmp-tree (:chr vc)
                                                       (:start-ci vc) (inc (:end-ci vc))))]
                (doseq [m-vc matches]
                  (reset! cmp-tree (remove-itree-vc @cmp-tree (:chr m-vc)
                                                    (:start m-vc) (:end m-vc))))
                [(if (seq matches) :sv-concordant (:1 disc-kwds)) (:vc vc)]))
            (remaining-cmp-svs [itree]
              (partition 2
                         (interleave (repeat (:2 disc-kwds)) (map :vc (get-itree-all itree)))))]

      (concat
       (map check-sv-concordance (parse-vcf-sv fname1 ref :interval-file interval-file
                                               :params params))
       (remaining-cmp-svs @cmp-tree)))))

(defn find-non-svs
  "Retrieve list of non-structural variants in the provided input file."
  [kwd vcf-source params]
  (->> (parse-vcf vcf-source)
       (filter #(nil? (get-sv-type % params)))
       (map :vc)
       (interleave (repeat kwd))
       (partition 2)))

(defn compare-sv
  "Compare structural variants, producing concordant and discordant outputs"
  [sample c1 c2 ref & {:keys [out-dir interval-file params]
                       :or {params {}}}]
  (let [base-out (str (fs/file (if (nil? out-dir) (fs/parent (:file c1)) out-dir)
                               (str sample "-%s-%s-%s.vcf")))
        disc-kwds {:1 (keyword (str "sv-" (:name c1) "-discordant"))
                   :2 (keyword (str "sv-" (:name c2) "-discordant"))}
        out-files (ordered-map
                   :sv-concordant (format base-out (:name c1) (:name c2) "svconcordance")
                   (:1 disc-kwds) (format base-out (:name c1) (:name c2) "svdiscordance")
                   (:2 disc-kwds) (format base-out (:name c2) (:name c1) "svdiscordance")
                   :nosv1 (itx/add-file-part (:file c1) "nosv" out-dir)
                   :nosv2 (itx/add-file-part (:file c2) "nosv" out-dir))]
    (when (itx/needs-run? (vals out-files))
      (with-open [vcf1-s (get-vcf-source (:file c1) ref :codec (structural-vcfcodec))
                  vcf2-s (get-vcf-source (:file c2) ref :codec (structural-vcfcodec))]
        (write-vcf-w-template (:file c1) out-files
                              (concat
                               (find-concordant-svs (:file c1) (:file c2) disc-kwds
                                                    ref interval-file params)
                               (find-non-svs :nosv1 vcf1-s params)
                               (find-non-svs :nosv2 vcf2-s params))
                              ref))
      ;; Remove SV VCF indexes since they use alternative Codecs
      (doseq [fname (vals out-files)]
        (let [x (str fname ".idx")]
          (if (fs/exists? x)
            (fs/delete x)))))
    out-files))

(defn compare-sv-pipeline
  "Handle input decomposition running SV detection through the standard pipeline."
  [c1 c2 exp config]
  (let [out-dir (get-in config [:dir :prep] (get-in config [:dir :out]))
        intervals (get c1 :intervals (get c2 :intervals (:intervals exp)))
        params (get exp :params {:min-indel 10 :default-cis [[100 10] [1000 200] [1e6 500]]})
        out-files (compare-sv (:sample exp) c1 c2 (:ref exp) :out-dir out-dir
                              :interval-file intervals :params params)]
    [(assoc c1 :file (:nosv1 out-files))
     (assoc c2 :file (:nosv2 out-files))
     (-> out-files
         (dissoc :nosv1)
         (dissoc :nosv2))]))

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
