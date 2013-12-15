(ns bcbio.variation.utils.cgmetrics
  "Add metrics from Complete Genomics masterVar file to a VCF.
  This updates a converted VCF from Complete Genomics with metrics information
  allowing assessment and filtering."
  (:import [org.broadinstitute.variant.variantcontext VariantContextBuilder]
           [org.broadinstitute.variant.vcf VCFHeader VCFInfoHeaderLine
            VCFHeaderLineCount VCFHeaderLineType])
  (:use [clojure.java.io]
        [ordered.set :only (ordered-set)]
        [bcbio.variation.normalize :only [hg19-map]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template
                                               get-vcf-iterator]])
  (:require [clojure.data.csv :as csv]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]))

(defn- get-masterVar-metrics
  "Get lookup dictionary of CG variant metrics by position."
  [in-file]
  (letfn [(variant-score [line name]
            (let [alleles ["1" "2"]]
              (/ (apply + (map #(Float/parseFloat (get line (format "allele%sVarScore%s" % name)))
                               alleles))
                 (count alleles))))
          (allele-balance [line]
            (/ (Float/parseFloat (get line "referenceAlleleReadCount"))
               (Float/parseFloat (get line "totalReadCount"))))]
    (with-open [rdr (reader in-file)]
      (let [csv-iter (drop-while #(< (count %) 3)
                                 (csv/read-csv rdr :separator \tab))
            header (first csv-iter)]
        (reduce (fn [coll xs]
                  (let [line (zipmap header xs)]
                    (assoc coll [(get hg19-map (get line "chromosome"))
                                 (inc (Integer/parseInt (get line "begin")))]
                           {:depth (get line "totalReadCount")
                            :qual-eaf (variant-score line "EAF")
                            :qual-vaf (variant-score line "VAF")
                            :ab (allele-balance line)})))
                {} (rest csv-iter))))))

(defn- add-cgmetrics-iter
  "Provide iterator of variants with CG metrics added"
  [vcf-source metrics]
  (letfn [(update-cgmetrics [vc x]
            (-> (VariantContextBuilder. (:vc vc))
                (.attributes (assoc (:attributes vc)
                               "DPCALL" (:depth x)
                               "AB" (:ab x)
                               "QUALEAF" (:qual-eaf x)
                               "QUALVAF" (:qual-vaf x)))
                .make))]
    (map (fn [vc]
           (if-let [cur-metrics (get metrics [(:chr vc) (:start vc)])]
             (update-cgmetrics vc cur-metrics)
             (:vc vc)))
         (parse-vcf vcf-source))))

(defn- add-cgmetrics-header
  "Add CG metrics definitions to the VCF input header."
  [_ header]
  (let [new #{(VCFInfoHeaderLine. "DPCALL" 1
                                  VCFHeaderLineType/Integer "Total depth used for calls")
              (VCFInfoHeaderLine. "QUALEAF" 1
                                  VCFHeaderLineType/Float
                                  "Variant quality under equal allele fraction model (EAF)")
              (VCFInfoHeaderLine. "QUALVAF" 1
                                  VCFHeaderLineType/Float
                                  "Variant quality under maximum likelihood variable allele fraction model (VAF)")
              (VCFInfoHeaderLine. "AB" 1
                                  VCFHeaderLineType/Float "Allele Balance")}]
    (VCFHeader. (apply ordered-set (concat (.getMetaDataInInputOrder header) new))
                (.getGenotypeSamples header))))

(defn add-cgmetrics
  "Add metrics from Complete Genomics masterVar file to VCF."
  [vcf-file mastervar-file ref-file & {:keys [out-dir]}]
  (let [out-file (fsp/add-file-part vcf-file "cgmetrics" out-dir)]
    (when (itx/needs-run? out-file)
      (with-open [vcf-iter (get-vcf-iterator vcf-file ref-file)]
        (write-vcf-w-template vcf-file {:out out-file}
                              (add-cgmetrics-iter vcf-iter
                                                  (get-masterVar-metrics mastervar-file))
                              ref-file :header-update-fn add-cgmetrics-header)))
    out-file))
