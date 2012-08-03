(ns bcbio.variation.utils.callsummary
  "Summarize a set of calls derived from multiple inputs to help with identifying filtering patterns."
  (:use [clojure.java.io]
        [bcbio.variation.config :only [load-config]]
        [bcbio.variation.variantcontext :only [parse-vcf get-vcf-iterator
                                               get-vcf-retriever variants-in-region]])
  (:require [clojure.string :as string]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]))

(defn- get-prepped-fname
  [call exp config]
  (let [dirname (get-in config [:dir :prep])]
    (str (file dirname (format "%s-%s-nomnp.vcf" (:sample exp) (:name call))))))

(defn report-vrn-summary
  "Report details on a variant based on items found in inputs."
  [wtr vc retriever fname-map]
  (let [hits (variants-in-region retriever (:chr vc) (:start vc) (:end vc))]
    (.write wtr (str (string/join ","
                                  [(:chr vc) (:start vc)
                                   (.getBaseString (:ref-allele vc))
                                   (string/join ";" (map #(.getBaseString %) (:alt-alleles vc)))
                                   (string/join ";" (sort (vec (set (map #(get fname-map (:fname %))
                                                                         hits)))))])
                     "\n"))))

(defn annotate-with-callsummary
  "Annotate input VCF with summary details from input files."
  [in-file config-file]
  (let [config (load-config config-file)
        exp (-> config :experiments first)
        orig-files (filter fs/exists? (map #(get-prepped-fname % exp config) (:calls exp)))
        fname-map (zipmap orig-files (map :name (:calls exp)))
        retriever (apply get-vcf-retriever (cons (:ref exp) orig-files))
        out-file (str (itx/file-root in-file) ".csv")]
    (with-open [vrn-iter (get-vcf-iterator in-file (:ref exp))
                wtr (writer out-file)]
      (doseq [vc (filter #(empty? (:filters %)) (parse-vcf vrn-iter))]
        (report-vrn-summary wtr vc retriever fname-map)))
    out-file))
