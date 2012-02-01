;; Support phased haplotype comparisons between variant calls.
;; This compares a phased set of calls versus haploid reference calls.
;;
;; The comparison logic is:
;; - Group calls into regions based on phasing
;; - For each phase region:
;;   - Determine which set of haploid alleles to compare with the reference
;;   - With each position in this haploid
;;      - Compare to reference allele
;;      - If mismatch and alternate allele matches reference, then phasing error
;;      - If mismatch and neither allele matches, then calling error

(ns bcbio.variation.phasing
  (:use [bcbio.variation.variantcontext :only [parse-vcf get-vcf-retriever]]))

(defn is-phased?
  "Check for phasing on a single genotype variant context."
  [vc]
  {:pre [(= 1 (count (:genotypes vc)))]}
  (-> vc :genotypes first :genotype .isPhased))

(defn parse-phased-haplotypes [in-vcf]
  "Separate phased haplotypes provided in diploid input genome."
  (lazy-seq
   (loop [vcs (parse-vcf in-vcf)
          cur-hap []
          all-haps []]
     ;; 3 conditions:
     ;; 1. Out of variants; add the current one to the list and done
     ;; 2. No current haplotype variants or phased with the previous variant:
     ;;    add to the current haplotype
     ;; 3. A new haplotype: add existing haplotype to list and create new
     (cond
      (nil? (first vcs)) (if (empty? cur-hap) all-haps (conj all-haps cur-hap))
      (or (empty? cur-hap)
          (is-phased? (first vcs))) (recur (rest vcs) (conj cur-hap (first vcs)) all-haps)
          :else (recur (rest vcs) [(first vcs)] (conj all-haps cur-hap))))))

(defn highest-count [xs]
  "Retrieve the item with the highest count in the supplied list.
  We break ties by sorting by the actual list items"
  (->> (frequencies xs)
       (sort-by val >)
       (partition-by second)
       first
       (sort-by first)
       ffirst))

(defn- get-alleles
  "Retrieve alleles for a single genotype variant context."
  [vc]
  {:pre [(= 1 (count (:genotypes vc)))]}
  (-> vc :genotypes first :alleles))

(defn- matching-allele
  "Determine allele index where the variant context matches haploid reference."
  [vc ref-vcs]
  {:pre [(every? #(= 1 (count (-> % :genotypes first :alleles))) ref-vcs)
         (= 1 (count (:genotypes vc)))]}
  (highest-count
   (map #(.indexOf (get-alleles vc) (-> % get-alleles first)) ref-vcs)))

(defn cmp-allele-to-ref
  "Compare the haploid allele of a variant against the reference call."
  [vc ref-vcs i]
  {:pre [(= 2 (count (get-alleles vc)))]}
  (let [ref-alleles (set (map #(-> % get-alleles first) ref-vcs))
        call-hap (nth (get-alleles vc) i)]
    (cond
     (empty? ref-alleles) :discordant
     (contains? ref-alleles call-hap) :concordant
     (some (partial contains? ref-alleles) (get-alleles vc)) :phasing-error
     :else :discordant)))

(defn get-variant-type [vcs]
  "Retrieve the type of a set of variants involved in a comparison.
  :indel -- insertions or deletions of more than 1bp
  :snp -- Single nucleotide changes or single basepair changes
  :unknown -- Other classs of variations (structural)"
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

(defn- nomatch-het-alt? [vc ref-vcs]
  "Determine if the variant has a non-matching heterozygous alternative allele."
  (let [match-allele-i (matching-allele vc ref-vcs)
        no-match-alleles (remove nil? (map-indexed
                                       (fn [i x] (if-not (= i match-allele-i) x))
                                       (get-alleles vc)))]
    (and (= "HET" (-> vc :genotypes first :type))
         (not-every? #(.isReference %) no-match-alleles))))

(defn- comparison-metrics [vc ref-vcs i]
  "Provide metrics for comparison of haploid allele to reference calls."
  {:comparison (cmp-allele-to-ref vc ref-vcs i)
   :variant-type (get-variant-type (cons vc ref-vcs))
   :nomatch-het-alt (nomatch-het-alt? vc ref-vcs)
   :vc (:vc vc)})

(defn- score-phased-region [vcs ref-fetch]
  "Provide scoring metrics for a phased region against a haplotype reference."
  (letfn [(ref-alleles [x]
            (ref-fetch (:chr x) (:start x) (:end x)))
          (ref-match-allele [x]
            (matching-allele x (ref-alleles x)))]
    (let [cmp-allele-i (highest-count (map ref-match-allele vcs))]
      (map #(comparison-metrics % (ref-alleles %) cmp-allele-i) vcs))))

(defn score-phased-calls [call-vcf ref-vcf]
  "Score a called VCF against reference based on phased regions."
  (let [ref-fetch (get-vcf-retriever ref-vcf)]
    (map #(score-phased-region % ref-fetch)
         (parse-phased-haplotypes call-vcf))))
