(ns bcbio.variation.haploid
  "Convert diploid variants from GATK into haploid calls based on genotype likelihoods.
  We assess diploid GATK calls based on the phred-normalized likelihood (PL). Lower variant
  PLs are likely to be true and included. The GATK documentation contains a detailed example
  of the format and interpretation:
  http://gatk.vanillaforums.com/discussion/36/understanding-unifiedgenotypers-vcf-files/p1"
  (:import [org.broadinstitute.sting.utils.variantcontext 
            VariantContextBuilder GenotypesContext GenotypeBuilder])
  (:use [clojure.java.io]
        [bcbio.variation.variantcontext :only [parse-vcf get-vcf-iterator write-vcf-w-template]])
  (:require [clojure.string :as string]
            [bcbio.run.itx :as itx]))

;; ## Convert diploid -> haploid

(defn- get-haploid-thresh
  "Threshold to include a heterozygous allele as a haploid homozygote variant.
  Based on type of variant: SNPs have lower threshold of inclusion."
  [vc]
  (case (:type vc)
    "SNP" 1e-5
    1e-50))

(defn get-likelihoods
  "Retrieve all likelihoods (PL) for genotype."
  [g & {:keys [no-convert]}]
  (when (.hasLikelihoods g)
    (let [in-map (-> (.getLikelihoods g) (.getAsMap (nil? no-convert)))]
      (zipmap (map #(.name %) (keys in-map)) (vals in-map)))))

(defn- get-haploid-genotype
  "Retrieve updated genotype with haploid allele."
  [vc]
  (let [g (-> vc :genotypes first)]
    (letfn [(maybe-variant-haploid [g vc]
              (when-let [variant-prob (get (get-likelihoods g) "HOM_VAR")]
                (when (> variant-prob (get-haploid-thresh vc))
                  (first (filter #(and (.isNonReference %) (.isCalled %))
                                   (.getAlleles g))))))
            (extract-mixed-allele [alleles]
              (let [ready (remove #(.isNoCall %) alleles)]
                (when (= 1 (count ready))
                  (first ready))))
            (get-haploid-allele [g vc]
              (case (:type g)
                "HOM_VAR" (first (:alleles g))
                "MIXED" (extract-mixed-allele (:alleles g))
                "HET" (maybe-variant-haploid (:genotype g) vc)
                nil))]
      (when-let [allele (get-haploid-allele g vc)]
        (doto (-> vc :vc .getGenotypes GenotypesContext/copy)
          (.replace (-> (GenotypeBuilder. (:genotype g))
                        (.alleles [allele])
                        .make)))))))

(defn- convert-to-haploid
  "Convert diploid allele to haploid variant."
  [vc]
  {:pre (= 1 (:num-samples vc))}
  (if-let [genotype (get-haploid-genotype vc)]
    [:haploid (-> (VariantContextBuilder. (:vc vc))
                  (.genotypes genotype)
                  (.make))]
    [:unchanged (:vc vc)]))

(defn diploid-calls-to-haploid
  "Convert set of diploid GATK calls on a haploid genome based on likelihoods."
  [vcf ref & {:keys [out-dir]}]
  (let [out-files {:haploid (itx/add-file-part vcf "haploid" out-dir)
                   :unchanged (itx/add-file-part vcf "nonhaploid" out-dir)}]
    (when (itx/needs-run? (vals out-files))
      (with-open [vcf-iter (get-vcf-iterator vcf ref)]
        (write-vcf-w-template vcf out-files
                              (map convert-to-haploid (parse-vcf vcf-iter))
                              ref)))
    (:haploid out-files)))

;; ## Examine diploid metrics

(defn write-het-variant-pls
  "Write phred likelihoods for het calls to be haploid variants."
  [vcf-file ref-file & attrs]
  (letfn [(get-pl [vc]
            (let [g (-> vc :genotypes first :genotype)]
              (when (.hasLikelihoods g)
                (let [in-map (-> (.getLikelihoods g) (.getAsMap true))]
                  (get (zipmap (map #(.name %) (keys in-map)) (vals in-map))
                       "HOM_VAR")))))]
  (let [out-file (str (itx/file-root vcf-file) "-het-pls.csv")]
    (with-open [vcf-iter (get-vcf-iterator vcf-file ref-file)
                wtr (writer out-file)]
      (doseq [val (->> (parse-vcf vcf-iter)
                       (filter #(= "HET" (-> % :genotypes first :type)))
                       (map get-pl)
                       (remove nil?))]
        (.write wtr (str (string/join "," (cons val attrs)) "\n"))))
    out-file)))

(defn -main
  [vcf ref]
  (diploid-calls-to-haploid vcf ref))
