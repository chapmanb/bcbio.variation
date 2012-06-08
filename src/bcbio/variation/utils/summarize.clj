(ns bcbio.variation.utils.summarize
  "Collapse a multi-sample VCF file into a CSV, R data.frame ready, parameter summary."
  (:use [clojure.java.io]
        [ordered.map :only [ordered-map]]
        [bcbio.variation.config :only [load-config]]
        [bcbio.variation.variantcontext :only [parse-vcf get-vcf-source]])
  (:require [clojure.data.csv :as csv]
            [incanter.stats :as istats]
            [bcbio.run.itx :as itx]))

(defn- flatten-vc-samples
  "Provide sample information from variant genotypes."
  [out vc attrs]
  (let [variant-types ["HET" "HOM_VAR"]]
    (letfn [(add-variant-totals [out gs]
              (let [counts (frequencies (map :type gs))]
                (reduce (fn [coll [k v]] (assoc coll k v))
                        out (map (fn [k] [k (get counts k 0)])
                                 variant-types))))
            (get-attr-avg [k gs]
              (istats/mean (->> gs
                                (filter #(contains? (set variant-types) (:type %)))
                                (map #(get-in % [:attributes k]))
                                (remove nil?)
                                (map #(Float/parseFloat %)))))
            (add-attr-avgs [out gs attrs]
              (reduce (fn [coll k] (assoc coll (str k "_sample_mean")
                                          (get-attr-avg k gs)))
                      out attrs))]
      (-> out
          (add-variant-totals (:genotypes vc))
          (add-attr-avgs (:genotypes vc) attrs)))))

(defn- flatten-vc-attrs
  "Extract attributes of interest from INFO field of variant."
  [out vc attrs]
  (reduce (fn [coll k] (assoc coll k (get-in vc [:attributes k])))
          out attrs))

(defn- flatten-vc
  "Provide tabular variant representation with provided attributes and sample information."
  [attrs sample-attrs vc]
  (-> (reduce (fn [coll k] (assoc coll k (get vc k)))
              (ordered-map) [:chr :start :id :type :qual])
      (flatten-vc-attrs vc attrs)
      (flatten-vc-samples vc sample-attrs)))

(defn vcf-to-table
  "Convert a VCF input to flattened CSV table with provided attributes."
  [vcf ref attrs sample-attrs]
  (let [out-file (str (itx/file-root vcf) "-variantsum.csv")]
    (when (itx/needs-run? out-file)
      (with-open [vcf-source (get-vcf-source vcf ref)
                  wtr (writer out-file)]
        (doseq [[i out] (map-indexed vector
                                     (map (partial flatten-vc attrs sample-attrs)
                                          (parse-vcf vcf-source)))]
          (when (= i 0)
            (csv/write-csv wtr [(map name (keys out))]))
          (csv/write-csv wtr [(vals out)]))))
    out-file))

(defn vcf-to-table-config
  "Prep a set of VCF to table conversions from input configuration file."
  [config-file]
  (let [config (load-config config-file)]
    (doall
     (flatten
      (for [exp (:experiments config)]
        (for [call (:calls exp)]
          (vcf-to-table (:file call) (:ref exp) (get-in call [:summary :attrs])
                        (get-in call [:summary :sample-attrs]))))))))

(defn -main [config-file]
  (vcf-to-table config-file))
