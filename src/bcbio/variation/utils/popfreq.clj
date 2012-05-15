(ns bcbio.variation.utils.popfreq
  "Associate population allele frequency with a list of variants.
  Annotates the original file with population frequencies based on rs IDs.
  Arguments:
    - original VCF file
    - attribute ID to include population frequencies in file (ie. GMAF)
    - Description of new frequency for VCF header
    - population VCF file
    - attribute ID to use for frequencies from population file (ie. AF)
    - reference genome FASTA file"
  (:import [org.broadinstitute.sting.utils.variantcontext VariantContextBuilder]
           [org.broadinstitute.sting.utils.codecs.vcf VCFHeader VCFInfoHeaderLine
            VCFHeaderLineCount VCFHeaderLineType])
  (:use [ordered.set :only [ordered-set]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template
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
  (println "Allele freqs" (count want-ids)) 
  (with-open [vcf-source (get-vcf-source vcf-file ref)]
    (reduce (fn [coll vc]
              (if (and (not (nil? (:id vc)))
                       (not (nil? (get (:attributes vc) freq-id)))
                       (contains? want-ids (:id vc)))
                (assoc coll (:id vc) (get (:attributes vc) freq-id))
                coll))
            {} (parse-vcf vcf-source))))

(defn add-pop-freqs
  "Lazy generator of variant contexts with added population frequencies."
  [vcf-source allele-freqs freq-id]
  (letfn [(update-allele-freq [vc new-freq]
            (-> (VariantContextBuilder. (:vc vc))
                (.attributes (assoc (:attributes vc) freq-id new-freq))
                .make))]
    (map #(update-allele-freq % (get allele-freqs (:id %) 0.0))
         (parse-vcf vcf-source))))

(defn- add-popfreq-header
  "Add new population frequency information to the VCF input header if needed."
  [freq-id freq-desc]
  (fn [header]
    (if (contains? (set (map #(when (= "INFO" (.getKey %))
                                (.getID %)) (.getMetaData header)))
                   freq-id) header
        (let [new #{(VCFInfoHeaderLine. freq-id 1
                                        VCFHeaderLineType/Float freq-desc)}]
          (VCFHeader. (apply ordered-set (concat (.getMetaData header) new))
                      (.getGenotypeSamples header))))))

(defn -main [orig-vcf freq-id freq-desc pop-vcf pop-freq-id ref-file]
  (let [out-file (itx/add-file-part orig-vcf "popfreq")
        allele-freqs (get-allele-freqs pop-vcf ref-file
                                       (get-rsids orig-vcf ref-file)
                                       pop-freq-id)]
    (with-open [vcf-source (get-vcf-source orig-vcf ref-file)]
      (write-vcf-w-template orig-vcf {:out out-file}
                            (add-pop-freqs vcf-source allele-freqs freq-id)
                            ref-file
                            :header-update-fn (add-popfreq-header freq-id freq-desc)))
    out-file))
