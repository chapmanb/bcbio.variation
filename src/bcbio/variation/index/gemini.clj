(ns bcbio.variation.index.gemini
  "Index and retrieve variant associated population genetic and disease data.
   Built on the Gemini framework: https://github.com/arq5x/gemini"
  (:use [ordered.map :only [ordered-map]]
        [bcbio.variation.web.db :only [get-sqlite-db]])
  (:require [clojure.java.jdbc :as sql]
            [clojure.java.shell :as shell]
            [clojure.string :as string]
            [bcbio.run.itx :as itx]))

(defn gemini-installed? []
  (let [info (try (shell/sh "gemini" "-h")
                  (catch java.io.IOException _
                    {:exit -1}))]
    (zero? (:exit info))))

(defn index-variant-file
  "Pre-index a variant file with gemini"
  [in-file ref-file & {:keys [re-index?]}]
  (when (gemini-installed?)
    (let [index-file (str (itx/file-root in-file) "-gemini.db")]
      (when (or (itx/needs-run? index-file) re-index?)
        (itx/with-tx-file [tx-index index-file]
          (shell/sh "gemini" "load" "-v" in-file tx-index)))
      index-file)))

(def ^{:doc "Gemini metrics to expose"}
  gemini-metrics
  (ordered-map "aaf_1kg_all" {:range [0.0 1.0]
                              :desc "1000 genomes allele frequency, all populations"}
               "gms_illumina" {:range [0.0 100.0]
                               :desc "Genome Mappability Score with an Illumina error model"}))

(defn available-metrics
  "Retrieve metrics available from Gemini."
  [_]
  (when gemini-installed?
    (map (fn [[k v]] (assoc v :id k)) gemini-metrics)))

(defn vc-attr-retriever
  "Retrieve metrics by name from a gemini index for provided VariantContexts."
  [in-file ref-file]
  (let [index-db (index-variant-file in-file ref-file)]
    (fn [vc attr]
      (when index-db
        (sql/with-connection (get-sqlite-db index-db)
          (sql/with-query-results rows
            [(str "SELECT " (name attr)
                  " FROM variants WHERE chrom = ? AND start = ? and ref = ?")
             (:chr vc) (:start vc) (.getBaseString (:ref-allele vc))]
            (get (first rows) (keyword (string/lower-case attr)))))))))

(defn get-raw-metrics
  "Retrieve table of Gemini metrics keyed on variant names."
  [in-file ref-file & {:keys [metrics]}]
  (when-let [index-db (index-variant-file in-file ref-file)]
    (let [plot-metrics (filter (partial contains? gemini-metrics)
                               (or metrics (map :id (available-metrics in-file))))]
      (sql/with-connection (get-sqlite-db index-db)
        (sql/with-query-results rows
          [(str "SELECT chrom, start, ref, " (string/join ", " plot-metrics)
                " FROM variants WHERE filter is NULL ORDER BY chrom, start")]
          (doall (map (fn [orig]
                        (reduce (fn [coll x]
                                  (assoc coll x (get orig (keyword (string/lower-case x)))))
                                {:id [(string/replace (:chrom orig) "chr" "") (inc (:start orig)) (:ref orig)]}
                                plot-metrics))
                      rows)))))))