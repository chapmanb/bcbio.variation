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
        [bcbio.variation.config :only [load-config]]
        [bcbio.variation.combine :only [combine-variants]]
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
  [vcf-file ref want-ids targets]
  (println "Retrieving allele freqs" (count want-ids) targets)
  (with-open [vcf-source (get-vcf-source vcf-file ref)]
    (reduce (fn [coll vc]
              (if (and (not (nil? (:id vc)))
                       (contains? want-ids (:id vc)))
                (assoc coll (:id vc) (zipmap (map :new-id targets)
                                             (map #(get-in vc [:attributes (:orig-id %)] 0.0)
                                                  targets)))
                coll))
            {} (parse-vcf vcf-source))))

(defn add-pop-freqs
  "Lazy generator of variant contexts with added population frequencies."
  [vcf-source allele-freqs ann-ids]
  (letfn [(update-allele-freq [vc new-freqs]
            (-> (VariantContextBuilder. (:vc vc))
                (.attributes (reduce (fn [coll cur-id]
                                       (assoc coll cur-id (get new-freqs cur-id 0.0)))
                                     (:attributes vc) ann-ids))
                .make))]
    (map #(update-allele-freq % (get allele-freqs (:id %) {}))
         (parse-vcf vcf-source))))

(defn- add-popfreq-header
  "Add new population frequency information to the VCF input header if needed."
  [new-ids]
  (letfn [(header-has-id? [header test-id]
            (contains? (set (map #(when (= "INFO" (.getKey %))
                                    (.getID %)) (.getMetaData header)))
                       test-id))]
    (fn [_ header]
      (let [new (->> new-ids
                     (remove #(header-has-id? header (:new-id %)))
                     (map #(VCFInfoHeaderLine. (:new-id %) 1
                                               VCFHeaderLineType/Float (:desc %)))
                     set)]
        (VCFHeader. (apply ordered-set (concat (.getMetaData header) new))
                    (.getGenotypeSamples header))))))

(defn- add-annotations
  [call ref out-dir]
  (let [orig-vcf (if (coll? (:file call))
                   (combine-variants (:file call) ref :out-dir out-dir :merge-type :full)
                   (:file call))
        out-file (itx/add-file-part orig-vcf "popfreq" out-dir)
        allele-freqs (get-allele-freqs (get-in call [:annotate :file]) ref
                                       (get-rsids orig-vcf ref)
                                       (get-in call [:annotate :targets]))]
    (if (itx/needs-run? out-file)
      (with-open [vcf-source (get-vcf-source orig-vcf ref)]
        (write-vcf-w-template orig-vcf {:out out-file}
                              (add-pop-freqs vcf-source allele-freqs
                                             (map :new-id (get-in call [:annotate :targets])))
                              ref
                              :header-update-fn
                              (add-popfreq-header (get-in call [:annotate :targets])))))
    out-file))

(defn -main [config-file]
  (let [config (load-config config-file)]
    (doseq [exp (:experiments config)]
      (doseq [call (:calls exp)]
        (add-annotations call (:ref exp) (get-in config [:dir :out]))))))
