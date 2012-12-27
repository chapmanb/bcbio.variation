(ns bcbio.variation.haploid
  "Convert diploid variants from GATK into haploid calls based on genotype likelihoods.
  We assess diploid GATK calls based on the phred-normalized likelihood (PL). Lower variant
  PLs are likely to be true and included. The GATK documentation contains a detailed example
  of the format and interpretation:
  http://gatkforums.broadinstitute.org/discussion/1268/how-should-i-interpret-vcf-files-produced-by-the-gatk"
  (:import [org.broadinstitute.sting.utils.variantcontext 
            VariantContextBuilder GenotypesContext GenotypeBuilder Allele])
  (:use [clojure.java.io]
        [bcbio.variation.variantcontext :only [parse-vcf get-vcf-iterator write-vcf-w-template]])
  (:require [clojure.string :as string]
            [bcbio.run.itx :as itx]))

;; ## Convert diploid -> haploid

(defn- get-haploid-thresh
  "Threshold to include a heterozygous allele as a haploid homozygote variant.
  Based on type of variant: SNPs have lower threshold of inclusion.
  Includes two thresholds: for possibly being homozygous variant and
  for not being homozygous reference."
  [vc]
  (case (:type vc)
    "SNP" {"HOM_VAR" 1e-5
           "HOM_REF" 1e-20}
    {"HOM_VAR" 1e-50
     "HOM_REF" 1e-200}))

(defn get-likelihoods
  "Retrieve all likelihoods (PL) for genotype."
  [g & {:keys [no-convert]}]
  (when (and (.hasLikelihoods g)
             (> (count (vec (.getAsVector (.getLikelihoods g)))) 1))
    (let [pl-vec (vec (.getAsVector (.getLikelihoods g)))]
      (if (= (count pl-vec) 2)
        (zipmap ["HOM_REF" "HOM_VAR"] pl-vec)
        (let [in-map (-> (.getLikelihoods g) (.getAsMap (nil? no-convert)))]
          (zipmap (map #(.name %) (keys in-map)) (vals in-map)))))))

(defn het-can-be-haploid?
  "Can we reasonably convert a heterozygote call to haploid based on likelihoods?
   We allow a het to pass when the homozygous variant prob is less than
   our defined thresholds, or our homozygous reference prob is greater."
  [g vc]
  (let [probs (get-likelihoods g)
        thresh (get-haploid-thresh vc)]
    (letfn [(passes-var? []
              (when-let [p (get probs "HOM_VAR")]
                (> p (get thresh "HOM_VAR"))))
            (passes-ref? []
              (when-let [p (get probs "HOM_REF")]
                (< p (get thresh "HOM_REF"))))]
      (or (passes-var?) (passes-ref?)))))

(defn- get-haploid-genotype
  "Retrieve updated genotype with haploid allele."
  [vc]
  (letfn [(maybe-variant-haploid [g vc]
            (when (het-can-be-haploid? g vc)
              (first (filter #(and (.isNonReference %) (.isCalled %))
                               (.getAlleles g)))))
          (extract-mixed-allele [alleles]
            (let [ready (remove #(.isNoCall %) alleles)]
              (when (= 1 (count ready))
                (first ready))))
          (get-haploid-allele [g vc]
            (case (:type g)
              "HOM_VAR" (first (:alleles g))
              "MIXED" (extract-mixed-allele (:alleles g))
              "HET" (maybe-variant-haploid (:genotype g) vc)
              nil))
          (add-haploid-genotype [context g]
            (let [allele (or (get-haploid-allele g vc) Allele/NO_CALL)]
              (doto context
                (.replace (-> (GenotypeBuilder. (:genotype g))
                              (.alleles [allele])
                              .make)))))]
    (reduce add-haploid-genotype
            (-> vc :vc .getGenotypes GenotypesContext/copy)
            (:genotypes vc))))

(defn- has-called-allele?
  "Check for at least one called Allele in a list of Genotypes."
  [genotypes]
  (not-every? #(.isNoCall %) (mapcat #(.getAlleles %) genotypes)))

(defn- convert-to-haploid
  "Convert diploid allele to haploid variant."
  [vc]
  (let [new-genotype (get-haploid-genotype vc)]
    (if (has-called-allele? new-genotype)
      [:haploid (-> (VariantContextBuilder. (:vc vc))
                    (.genotypes new-genotype)
                    (.make))]
      [:unchanged (:vc vc)])))

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
