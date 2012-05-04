(ns bcbio.variation.phasing
  "Support phased haplotype comparisons between variant calls.
   Compares a phased set of calls versus haploid reference calls.

   The comparison logic is:

   - Group calls into regions based on phasing
   - For each phase region:
     - Determine which set of haploid alleles to compare with the reference
     - With each position in this haploid:
       - Compare to reference allele
       - If mismatch and alternate allele matches reference, then phasing error
       - If mismatch and neither allele matches, then calling error"
  (:import [org.broadinstitute.sting.utils.interval IntervalUtils IntervalSetRule]
           [org.broadinstitute.sting.utils GenomeLocParser GenomeLoc])
  (:use [bcbio.variation.variantcontext :only [parse-vcf get-vcf-retriever get-vcf-source
                                               write-vcf-w-template]]
        [bcbio.align.ref :only [get-seq-dict]]
        [bcbio.variation.callable :only [get-bed-source]]
        [ordered.map :only [ordered-map]])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]))

;; ## Find phased haplotypes in VCF

(defn is-phased?
  "Check for phasing on a single genotype variant context."
  [vc]
  {:pre [(= 1 (count (:genotypes vc)))]}
  (let [g (-> vc :genotypes first)]
    (or (= 1 (count (:alleles g)))
        (.isPhased (:genotype g)))))

(defn parse-phased-haplotypes
  "Separate phased haplotypes provided in diploid input genome.
   3 conditions:

   1. Out of variants; add the current one to the list and done
   2. No current haplotype variants or phased with the previous variant:
      add to the current haplotype
   3. A new haplotype: add existing haplotype to list and create new"
  [vcf-source]
  (lazy-seq
   (loop [vcs (parse-vcf vcf-source)
          cur-hap []
          all-haps []]
     (cond
      (nil? (first vcs)) (if (empty? cur-hap) all-haps (conj all-haps cur-hap))
      (or (empty? cur-hap)
          (is-phased? (first vcs))) (recur (rest vcs) (conj cur-hap (first vcs)) all-haps)
          :else (recur (rest vcs) [(first vcs)] (conj all-haps cur-hap))))))

;; ## Compare phased variants

(defn highest-count
  "Retrieve the item with the highest count in the supplied list.
  We break ties by sorting by the actual list items"
  [xs]
  (->> (frequencies xs)
       (sort-by val >)
       (partition-by second)
       first
       (sort-by first)
       ffirst))

(defn- get-alleles
  "Convenience function to get alleles for a single genotype variant context."
  [vc]
  {:pre [(= 1 (count (:genotypes vc)))]}
  (-> vc :genotypes first :alleles))

(defn- matching-allele
  "Determine allele index where the variant context matches haploid reference."
  [vc ref-vcs]
  {:pre [(every? #(= 1 (count (-> % :genotypes first :alleles))) ref-vcs)
         (= 1 (count (:genotypes vc)))]}
  (if (empty? ref-vcs)
    (.indexOf (get-alleles vc) (:ref-allele vc))
    (highest-count
     (remove neg?
             (map #(.indexOf (get-alleles vc) (-> % get-alleles first)) ref-vcs)))))

(defn cmp-allele-to-ref
  "Compare the haploid allele of a variant against the reference call."
  [vc ref-vcs i]
  (letfn [(is-ref-allele? [x]
            (= (.getBaseString x) (-> vc :ref-allele .getBaseString)))]
    (let [ref-alleles (set (map #(-> % get-alleles first) ref-vcs))
          call-hap (when-not (nil? i) (nth (get-alleles vc) i))]
      (cond
       (nil? call-hap) :discordant
       (and (is-ref-allele? call-hap)
            (or (empty? ref-alleles)
                (contains? ref-alleles call-hap))) :ref-concordant
       (empty? ref-alleles) :discordant
       (contains? ref-alleles call-hap) :concordant
       (some (partial contains? ref-alleles) (get-alleles vc)) :phasing-error
       :else :discordant))))

(defn get-variant-type
  "Retrieve the type of a set of variants involved in a comparison.

  - `:indel` -- insertions or deletions of more than 1bp
  - `:snp` -- Single nucleotide changes or single basepair changes
  - `:unknown` -- Other classs of variations (structural)"
  [vcs]
  (letfn [(is-indel? [x]
            (= "INDEL" (:type x)))
          (is-multi-indel? [x]
            (and (is-indel? x)
                 (not-every? #(contains? #{0 1} %)
                             (map #(-> % .getBaseString count) (get-alleles x)))))
          (is-snp? [x]
            (= "SNP" (:type x)))]
    (cond
     (some is-multi-indel? vcs) :indel
     (some is-indel? vcs) :snp
     (every? is-snp? vcs) :snp
     :else :unknown)))

(defn- nomatch-het-alt?
  "Determine if the variant has a non-matching heterozygous alternative allele."
  [vc ref-vcs]
  (let [match-allele-i (matching-allele vc ref-vcs)
        no-match-alleles (remove nil? (map-indexed
                                       (fn [i x] (if-not (= i match-allele-i) x))
                                       (get-alleles vc)))]
    (and (= "HET" (-> vc :genotypes first :type))
         (not-every? #(.isReference %) no-match-alleles))))

(defn- comparison-metrics
  "Provide metrics for comparison of haploid allele to reference calls."
  [vc ref-vcs i]
  {:comparison (cmp-allele-to-ref vc ref-vcs i)
   :variant-type (get-variant-type (cons vc ref-vcs))
   :nomatch-het-alt (nomatch-het-alt? vc ref-vcs)
   :vc (:vc vc)
   :ref-vcs (map :vc ref-vcs)})

(defn- score-phased-region
  "Provide scoring metrics for a phased region against a haplotype reference."
  [vcs ref-fetch]
  (letfn [(get-ref-vcs [x]
            (ref-fetch (:chr x) (:start x) (:end x)))
          (ref-match-allele [x]
            (matching-allele x (get-ref-vcs x)))]
    (let [cmp-allele-i (highest-count (map ref-match-allele vcs))]
      (map #(comparison-metrics % (get-ref-vcs %) cmp-allele-i) vcs))))

(defn score-phased-calls
  "Score a called VCF against reference based on phased regions."
  [call-vcf-s ref-vcf-s]
  (let [ref-fetch (get-vcf-retriever ref-vcf-s)]
    (map #(score-phased-region % ref-fetch)
         (parse-phased-haplotypes call-vcf-s))))

;; ## Summarize phased comparisons

(defn- write-concordance-output
  "Write concordant and discordant variants to VCF output files."
  [vc-info to-capture sample-name base-info other-info out-dir ref]
  (let [base-dir (if (nil? out-dir) (fs/parent (:file base-info)) out-dir)
        gen-file-name (fn [x] (str (fs/file base-dir (format "%s-%s-%s-%s.vcf"
                                                             sample-name (:name base-info)
                                                             (:name other-info) (name x)))))
        out-files (apply ordered-map (flatten (map (juxt identity gen-file-name)
                                                   to-capture)))]
    (if-not (fs/exists? base-dir)
      (fs/mkdirs base-dir))
    (when (itx/needs-run? (vals out-files))
      (write-vcf-w-template (:file base-info) out-files
                            (filter #(contains? (set to-capture) (first %))
                                    (map (juxt :comparison :vc)
                                         (flatten vc-info)))
                            ref))
    out-files))

(defn count-comparison-bases
  "Provide counts for comparison: entire region plus user specified regions"
  [total-bed call-bed ref-file]
  (letfn [(feature-size [x]
            (cond
             (instance? GenomeLoc x) (- (.getStop x) (.getStart x))
             :else (- (.getEnd x) (.getStart x))))
          (count-bases [xs]
            (apply + (map feature-size xs)))
          (genome-loc-list [x]
            (let [parser (GenomeLocParser. (get-seq-dict ref-file))]
              (with-open [bed-source (get-bed-source x)]
                (->> bed-source
                     .iterator
                     (map #(.createGenomeLoc parser %))
                     doall))))
          (merge-intervals [x y]
            (IntervalUtils/mergeListsBySetOperator (genome-loc-list x)
                                                   (genome-loc-list y)
                                                   IntervalSetRule/INTERSECTION))]
    (with-open [bed-source (get-bed-source total-bed)]
      (let [total (count-bases (.iterator bed-source))
            compared (if (nil? call-bed) total
                         (count-bases (merge-intervals total-bed call-bed)))]
        {:percent (* 100.0 (/ compared total))
         :compared compared
         :total total}))))

(defn- get-phasing-metrics
  "Collect summary metrics for concordant/discordant and phasing calls"
  [vc-info exp-interval-file call-interval-file ref-file]
  (letfn [(count-nomatch-het-alt [xs]
            (count (filter #(and (contains? #{:concordant :ref-concordant} (:comparison %))
                                 (:nomatch-het-alt %))
                           (flatten vc-info))))
          (blank-count-dict []
            {:snp 0 :indel 0})
          (add-current-count [coll x]
            (let [cur-val (map x [:comparison :variant-type])]
              (assoc-in coll cur-val (inc (get-in coll cur-val)))))]
    (reduce add-current-count
            {:haplotype-blocks (count vc-info)
             :total-bases (count-comparison-bases exp-interval-file call-interval-file ref-file)
             :nonmatch-het-alt (count-nomatch-het-alt vc-info)
             :concordant (blank-count-dict)
             :ref-concordant (blank-count-dict)
             :discordant (blank-count-dict)
             :phasing-error (blank-count-dict)}
            (flatten vc-info))))

;; ## Entry point for phased haploid VCF comparisons

(defmulti compare-two-vcf-phased
  "Compare two VCF files including phasing with a haplotype reference
  Handle grading special case as well as standard comparisons."
  (fn [_ exp _] (keyword (get exp :approach "compare"))))

(defmethod compare-two-vcf-phased :grade
  [phased-calls exp config]
  {:pre [(= 1 (count (get phased-calls true)))
         (= 1 (count (get phased-calls false)))]}
  (let [ref (first (get phased-calls true))
        call (first (get phased-calls false))]
    (with-open [ref-vcf-s (get-vcf-source (:file ref) (:ref exp))
                call-vcf-s (get-vcf-source (:file call) (:ref exp))]
      (let [compared-calls (score-phased-calls call-vcf-s ref-vcf-s)]
        {:c-files (write-concordance-output compared-calls
                                            [:concordant :discordant :phasing-error]
                                            (:sample exp) call ref
                                            (get-in config [:dir :out]) (:ref exp))
         :metrics (get-phasing-metrics compared-calls (:intervals exp)
                                       (:intervals call) (:ref exp))
         :c1 call :c2 ref :sample (:sample exp) :exp exp}))))

(defn- convert-vcs-to-compare
  "Convert stream of variant context haploid comparison to standard,
  keyed by :concordant and :discordant-name keywords."
  [name1 name2 cmps]
  (letfn [(update-keyword [x]
            (let [new-xs (remove #(< (.getStart (:vc %)) (.getStart (:vc x)))
                                 (sort-by #(get-in % [:vc :start])
                                          (map #(-> x (assoc :vc %) (dissoc :ref-vcs))
                                               (:ref-vcs x))))
                  [dis-kw1 dis-kw2] (map #(keyword (format "%s-discordant" %)) [name1 name2])]
              (case (:comparison x)
                :concordant new-xs
                (:discordant :phasing-error) (cons
                                              (assoc x :comparison dis-kw2)
                                              (map #(assoc % :comparison dis-kw1) new-xs))
                nil)))]
    (remove nil?
            (flatten
             (map update-keyword (flatten cmps))))))

(defmethod compare-two-vcf-phased :compare
  [phased-calls exp config]
  {:pre [(= 2 (count (flatten (vals phased-calls))))
         (pos? (count (get phased-calls true)))]}
  (let [cmp1 (first (get phased-calls true))
        cmp2 (if-let [nophased (get phased-calls false)]
               (first nophased)
               (second (get phased-calls true)))
        to-capture (concat [:concordant]
                           (map #(keyword (format "%s-discordant" (:name %)))
                                [cmp1 cmp2]))]
    (with-open [vcf1-s (get-vcf-source (:file cmp1) (:ref exp))
                vcf2-s (get-vcf-source (:file cmp2) (:ref exp))]
      (let [compared-calls (convert-vcs-to-compare (:name cmp1) (:name cmp2)
                                                   (score-phased-calls vcf2-s vcf1-s))]
        {:c-files (write-concordance-output compared-calls to-capture
                                            (:sample exp) cmp1 cmp2 
                                            (get-in config [:dir :out]) (:ref exp))
         :c1 cmp1 :c2 cmp2 :sample (:sample exp) :exp exp}))))

;; ## Utility functions

(defn is-haploid?
  "Is the provided VCF file a haploid genome (one genotype or all homozygous).
  Samples the first set of variants, checking for haploid calls."
  [vcf-file ref-file]
  (let [sample-size 1000]
    (letfn [(is-vc-haploid? [vc]
              (or (= 1 (apply max (map #(count (:alleles %)) (:genotypes vc))))
                  (contains? #{"HOM_REF" "HOM_VAR"} (:type vc))))]
      (with-open [vcf-source (get-vcf-source vcf-file ref-file)]
        (every? is-vc-haploid? (take sample-size (parse-vcf vcf-source)))))))
