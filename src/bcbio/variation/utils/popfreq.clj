(ns bcbio.variation.utils.popfreq
  "Associate population allele frequency with a list of variants.
  Annotates the original file with population frequencies if found."
  (:import [org.broadinstitute.sting.utils.variantcontext VariantContextBuilder])
  (:use [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template
                                               get-vcf-source]])
  (:require [bcbio.run.itx :as itx]))

(defn get-rsids
  "Retrieve all rsIDs from the input vcf-file"
  [vcf-file ref]
  (with-open [vcf-source (get-vcf-source vcf-file ref)]
    (set (remove nil? (map :id (parse-vcf vcf-source))))))

(defn get-allele-freqs
  "Retrieve allele frequencies from population VCF for IDs of interest."
  [vcf-file ref want-ids freq-id]
  (with-open [vcf-source (get-vcf-source vcf-file ref)]
    (reduce (fn [coll vc]
              (if (and (contains? want-ids (:id vc))
                       (not (nil? (get (:attributes vc) freq-id))))
                (assoc coll (:id vc) (get (:attributes vc) freq-id))
                coll))
            {} (parse-vcf vcf-source))))

(defn add-pop-freqs
  "Lazy generator of variant contexts with added population frequencies."
  [vcf-source allele-freqs]
  (letfn [(update-allele-freq [vc new-freq]
            (-> (VariantContextBuilder. (:vc vc))
                (.attributes (assoc (:attributes vc) "AF" new-freq))
                .make))]
    (remove nil?
            (map #(when-let [freq (get allele-freqs (:id %))]
                    (update-allele-freq % freq))
                 (parse-vcf vcf-source)))))

(defn -main [orig-vcf pop-vcf ref-file]
  (let [freq-id "AF"
        out-file (itx/add-file-part orig-vcf "popfreq")
        allele-freqs (get-allele-freqs pop-vcf ref-file
                                       (get-rsids orig-vcf ref-file)
                                       freq-id)]
    (with-open [vcf-source (get-vcf-source orig-vcf ref-file)]
      (write-vcf-w-template orig-vcf {:out out-file}
                            (add-pop-freqs vcf-source allele-freqs)
                            ref-file))
    out-file))
