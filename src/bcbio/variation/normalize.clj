(ns bcbio.variation.normalize
  "Prepare a VCF file for comparison by normalizing chromosome names,
  sort order, sample name, and genotype representation.
  This handles the work of making slightly different representations
  match, enabling VCF comparisons."
  (:use [bcbio.variation.variantcontext :only (parse-vcf write-vcf-w-template
                                               get-seq-dict vcf-source)]))

(defn prep-vcf
  "Prepare VCF for comparison by normalizing high level attributes"
  [in-vcf ref-file sample]
  (let [ref-chrs (map #(.getSequenceName %) (-> ref-file get-seq-dict .getSequences))
        vcf (vcf-source in-vcf)
        vcf-chrs (.getSequenceNames vcf)]
    (println ref-chrs)
    (println vcf-chrs)))
