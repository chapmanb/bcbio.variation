(ns bcbio.variation.index.metrics
  "Pre-index a variant file for quick retrieval of associated metrics."
  (:use [ordered.map :only [ordered-map]]
        [bcbio.variation.filter.classify :only [get-vc-attrs]]
        [bcbio.variation.variantcontext :only [get-vcf-header get-vcf-iterator parse-vcf]]
        [bcbio.variation.web.db :only [get-sqlite-db]])
  (:require [clojure.string :as string]
            [clojure.java.jdbc :as sql]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]))

(def ^{:doc "Metrics to expose, ranked in order of priority with default min/max values."}
  expose-metrics
  (ordered-map "QUAL" [0.0 100000.0]
               "DP" [0.0 5000.0]
               "MQ" [0.0 75.0]
               "QD" [0.0 200.0]
               "HaplotypeScore" [0.0 250.0]))

(def ^{:doc "Default metrics that are always available." :private true}
  default-metrics
  [{:id "QUAL" :desc "Variant quality score, phred-scaled"}])

(defn available-metrics
  "Retrieve metrics available for variant input file."
  [vcf-file]
  (letfn [(convert-header [line]
            {:id (.getID line)
             :desc (.getDescription line)})]
    (let [metrics-order (reduce (fn [coll [i x]] (assoc coll x i))
                                {} (map-indexed vector (keys expose-metrics)))]
      (->> (get-vcf-header vcf-file)
           .getMetaDataInInputOrder
           (filter #(= "INFO" (.getKey %)))
           (filter #(contains? expose-metrics (.getID %)))
           (map convert-header)
           (concat default-metrics)
           (sort-by #(get metrics-order (:id %)))))))

(defn- create-metrics-tables
  "Create table to represent variant metrics"
  [metrics]
  (apply sql/create-table (concat [:metrics
                                   [:contig :text]
                                   [:start :integer]
                                   [:refallele :text]]
                                  (map (fn [x] [(:id x) :float]) metrics))))

(defn index-variant-file
  "Pre-index a variant file with associated metrics."
  [in-file ref-file]
  (let [metrics (available-metrics in-file)
        index-file (str (itx/file-root in-file) "-metrics.db")]
    (when-not (fs/exists? index-file)
      (sql/with-connection (get-sqlite-db index-file :create true)
        (sql/transaction
         (create-metrics-tables metrics))
        (with-open [vcf-iter (get-vcf-iterator in-file ref-file)]
          (doseq [vc (parse-vcf vcf-iter)]
            (sql/transaction
             (sql/insert-record :metrics
                                (-> (get-vc-attrs vc (map :id metrics))
                                    (assoc :contig (:chr vc))
                                    (assoc :start (:start vc))
                                    (assoc :refallele (.getBaseString (:ref-allele vc))))))))))
    index-file))

(defn get-raw-metrics
  "Retrieve table of raw metrics using indexed variant file"
  [in-file ref-file & {:keys [metrics]}]
  (let [index-db (index-variant-file in-file ref-file)
        plot-metrics (if (nil? metrics) (map :id (available-metrics in-file)) metrics)]
    (sql/with-connection (get-sqlite-db index-db)
      (sql/with-query-results rows
        [(str "SELECT contig, start, refallele, " (string/join ", " plot-metrics)
              " FROM metrics ORDER BY contig, start")]
        (doall (map (fn [orig]
                      (reduce (fn [coll x]
                                (assoc coll x (get orig (keyword (string/lower-case x)))))
                              {:id (format "%s-%s-%s" (:contig orig) (:start orig)
                                           (:refallele orig))}
                              plot-metrics))
                    rows))))))