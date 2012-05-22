(ns bcbio.variation.haploid
  "Convert diploid variants from GATK into haploid calls based on genotype likelihoods.
  We assess diploid GATK calls based on the phred-normalized likelihood (PL). Lower variant
  PLs are likely to be true and included. The GATK documentation contains a detailed example
  of the format and interpretation:
  http://www.broadinstitute.org/gsa/wiki/index.php/
  Understanding_the_Unified_Genotyper%27s_VCF_files#How_genotypes_are_represented_in_a_VCF"
  (:import [org.broadinstitute.sting.utils.variantcontext 
            VariantContextBuilder GenotypesContext Genotype])
  (:use [bcbio.variation.variantcontext :only [parse-vcf get-vcf-source write-vcf-w-template]])
  (:require [bcbio.run.itx :as itx]))

(def ^{:doc "Threshold to include a heterozygous allele as a haploid homozygote variant."}
  haploid-thresh 1e-5)

(defn- get-haploid-genotype
  "Retrieve updated genotype with haploid allele."
  [vc]
  (let [g (-> vc :genotypes first)]
    (letfn [(maybe-variant-haploid [g]
              (when (.hasLikelihoods g)
                (let [in-map (-> (.getLikelihoods g) (.getAsMap true))
                      variant-prob (get (zipmap (map #(.name %) (keys in-map)) (vals in-map))
                                        "HOM_VAR")]
                  (when (> variant-prob haploid-thresh)
                    (first (filter #(and (.isNonReference %) (.isCalled %))
                                   (.getAlleles g)))))))
            (extract-mixed-allele [alleles]
              (let [ready (remove #(.isNoCall %) alleles)]
                (when (= 1 (count ready))
                  (first ready))))
            (get-haploid-allele [g]
              (case (:type g)
                "HOM_VAR" (first (:alleles g))
                "MIXED" (extract-mixed-allele (:alleles g))
                "HET" (maybe-variant-haploid (:genotype g))
                nil))]
      (when-let [allele (get-haploid-allele g)]
        (doto (-> vc :vc .getGenotypes GenotypesContext/copy)
          (.replace (Genotype/modifyAlleles (:genotype g) [allele])))))))

(defn- convert-to-haploid
  "Convert diploid allele to haploid variant."
  [vc]
  {:pre (= 1 (count (:genotypes vc)))}
  (when-let [genotype (get-haploid-genotype vc)]
    (-> (VariantContextBuilder. (:vc vc))
        (.genotypes genotype)
        (.make))))

(defn diploid-calls-to-haploid
  "Convert set of diploid GATK calls on a haploid genome based on likelihoods."
  [vcf ref & {:keys [out-dir]}]
  (let [out-files {:out (itx/add-file-part vcf "haploid" out-dir)}]
    (when (itx/needs-run? (vals out-files))
      (with-open [vcf-source (get-vcf-source vcf ref)]
        (write-vcf-w-template vcf out-files
                              (remove nil?
                                      (map convert-to-haploid (parse-vcf vcf-source)))
                              ref)))
    (:out out-files)))
