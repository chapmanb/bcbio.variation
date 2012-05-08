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
  (:use [bcbio.variation.structural :only [prep-itree get-itree-overlap
                                           remove-itree-vc get-itree-all]]
        [bcbio.variation.variantcontext :only [parse-vcf get-vcf-retriever get-vcf-source
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

(defn cmp-allele-to-expected
  "Compare the haploid allele of a variant against the expected call."
  [vc e-vc i]
  (letfn [(is-ref-allele? [x]
            (apply = (map #(.getBaseString (% x)) [:cmp :ref])))
          (get-cmp-allele [i x]
            {:ref (:ref-allele x)
             :cmp (nth (get-alleles x) i)})
          (get-all-alleles [x]
            (map #(get-cmp-allele % x) (range (count (get-alleles x)))))]
    (let [e-allele (when-not (nil? e-vc)
                       (get-cmp-allele 0 e-vc))
          call-hap (when-not (or (nil? i) (nil? vc) (neg? i))
                     (get-cmp-allele i vc))]
      (cond
       (nil? call-hap) :discordant
       (and (is-ref-allele? call-hap)
            (or (nil? e-allele)
                (= e-allele call-hap))) :ref-concordant
       (nil? e-allele) :discordant
       (= e-allele call-hap) :concordant
       (some (partial = e-allele) (get-all-alleles vc)) :phasing-error
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
  [vc e-vc]
  (let [match-allele-i (matching-allele vc [e-vc])
        no-match-alleles (remove nil? (map-indexed
                                       (fn [i x] (if-not (= i match-allele-i) x))
                                       (get-alleles vc)))]
    (and
         (= "HET" (-> vc :genotypes first :type))
         (not-every? #(.isReference %) no-match-alleles))))

(defn- comparison-metrics
  "Provide metrics for comparison of haploid expected alleles to variant calls."
  [cmp-itree i e-vc]
  (let [cmp-vc (->> (get-itree-overlap cmp-itree (:chr e-vc) (:start e-vc) (:end e-vc))
                    (filter #(= (:start %) (:start e-vc)))
                    first)]
    {:comparison (cmp-allele-to-expected cmp-vc e-vc i)
     :variant-type (get-variant-type [cmp-vc e-vc])
     :nomatch-het-alt (when (not-any? nil? [cmp-vc e-vc]) (nomatch-het-alt? cmp-vc e-vc))
     :start (if (nil? cmp-vc) (:start e-vc) (:start cmp-vc))
     :end (when-not (nil? cmp-vc) (:end cmp-vc))
     :vc (:vc cmp-vc)
     :ref-vc (:vc e-vc)}))

(defn- score-phased-region
  "Provide scoring metrics for a phased region against expected haplotype variants."
  [expect-fetch vcs]
  (letfn [(get-ref-vcs [x]
            (expect-fetch (:chr x) (:start x) (:end x)))
          (ref-match-allele [x]
            (matching-allele x (expect-fetch (:chr x) (:start x) (:end x))))
          (get-regional-expected-vcs
            [itree]
            {:pre [(= 1 (count (keys itree)))]}
            (let [[chr tree] (first itree)]
              (sort-by :start
                       (expect-fetch chr (-> tree .min .getStart)
                                     (dec (-> tree .max .getEnd))))))
          (compare-and-update-itree [info e-vc]
            (let [cmp (comparison-metrics (:itree info) (:cmp-i info) e-vc)]
              (-> info
                  (assoc :itree (remove-itree-vc (:itree info) (:chr e-vc)
                                                 (:start cmp) (:end cmp)))
                  (assoc :out (cons cmp (:out info))))))
          (add-unmapped-cmps [info]
            (concat (:out info)
                    (map (fn [vc] {:comparison (cmp-allele-to-expected vc nil (:cmp-i info))
                                   :variant-type (get-variant-type [vc])
                                   :nomatch-het-alt false
                                   :start (:start vc)
                                   :vc (:vc vc)
                                   :ref-vc nil})
                         (get-itree-all (:itree info)))))]
    (let [cmp-allele-i (highest-count (map ref-match-allele vcs))
          vc-itree (prep-itree vcs :start :end)]
      (->> (reduce compare-and-update-itree
                   {:itree vc-itree :cmp-i cmp-allele-i :out []}
                   (get-regional-expected-vcs vc-itree))
           add-unmapped-cmps
           (sort-by :start)))))

(defn score-phased-calls
  "Score a called VCF against expected haploid variants based on phased regions."
  [call-vcf-s ref-vcf-s]
  (map (partial score-phased-region (get-vcf-retriever ref-vcf-s))
       (parse-phased-haplotypes call-vcf-s)))

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
            (let [ref-x (-> x
                            (assoc :vc (:ref-vc x))
                            (dissoc :ref-vc))
                  [dis-kw1 dis-kw2] (map #(keyword (format "%s-discordant" %)) [name1 name2])]
              (case (:comparison x)
                :concordant [ref-x]
                (:discordant :phasing-error) [(assoc x :comparison dis-kw2)
                                              (assoc ref-x :comparison dis-kw1)]
                nil)))]
    (remove #(or (nil? %) (nil? (:vc %)))
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
