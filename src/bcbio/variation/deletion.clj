(ns bcbio.variation.deletion
  "Provide custom comparisons for deletions in reference genomes.
  Deletion comparisons are more difficult since they remove additional
  sequence for future comparison."
  (:use [bcbio.variation.variantcontext :only [parse-vcf get-vcf-retriever get-vcf-source]]))

(defn is-deletion?
  [vc]
  (and (= (:type vc) "INDEL")
       (pos? (.length (:ref-allele vc)))))

(defn compare-deletion-vc
  "Compare a deletion at a variant context, returning match information."
  [cmp-fetch vc]
  (let [cmps (filter is-deletion? (apply cmp-fetch (map #(% vc) [:chr :start :end])))]
    {:vc vc
     :ref-vcs cmps
     :comparison (cond
                  (empty? cmps) :discordant)}))

(defn compare-deletions
  "Compare deletions in two VCF files, returning concordant and discordant calls."
  [ref-vcf cmp-vcf ref]
  (with-open [ref-vcf-s (get-vcf-source ref-vcf ref)
              cmp-vcf-s (get-vcf-source cmp-vcf ref)]
    (let [cmp-fetch (get-vcf-retriever cmp-vcf-s)]
      (doall
       (map (partial compare-deletion-vc cmp-fetch)
            (filter is-deletion? (parse-vcf ref-vcf-s)))))))
