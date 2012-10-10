(ns bcbio.variation.api.metrics
  "Provide high level API for accessing variant associated metrics."
  (:import [org.jfree.data.statistics HistogramDataset HistogramType])
  (:use [bcbio.variation.api.shared :only [web-config]]
        [bcbio.variation.variantcontext :only [get-vcf-iterator parse-vcf]])
  (:require [bcbio.variation.index.metrics :as im]
            [bcbio.variation.index.gemini :as gemini]
            [bcbio.variation.remote.core :as remote]))

;; ## Helper functions

(declare available-metrics)

(defn- get-histogram-bins
  [items n bin-min bin-max]
  "Retrieve values binned into a histogram using JFree Chart."
  (let [ds (doto (HistogramDataset.)
             (.setType HistogramType/RELATIVE_FREQUENCY)
             (.addSeries 0 (double-array items) n bin-min bin-max))]
    {:x (map #(.getXValue ds 0 %) (range (.getItemCount ds 0)))
     :y (map #(.getYValue ds 0 %) (range (.getItemCount ds 0)))}))

(defn- clean-raw-metrics
  "Remove nil values and empty input metrics."
  [raw metrics]
  (reduce (fn [coll [k vs]]
            (let [clean-vs (remove nil? vs)]
              (if (empty? clean-vs)
                coll
                (assoc coll k clean-vs))))
          {}
          (zipmap metrics (map (fn [x] (map #(get % x) raw)) metrics))))

(defn- get-metric-range
  "Retrieve the configured range of a metric"
  [metric]
  (-> (filter #(= (:id %) metric) (available-metrics nil))
      first
      :range))

(defn- prepare-plot-metrics
  "Bin metrics in preparation for histogram display using predefined min-max boundaries."
  [metric raw]
  (let [bins 100
        [bin-min bin-max] (get-metric-range metric)
        data (get-histogram-bins raw bins bin-min bin-max)]
    {:vals (:y data)
     :bin-width (- (second (:x data)) (first (:x data)))
     :x-scale {:type :linear
               :domain [bin-min bin-max]}
     :y-scale {:type :linear}}))

(defn- combined-raw-metrics
  "Retrieve raw metrics from multiple sources, combining on IDs."
  [vcf-file ref-file metrics use-subsample?]
  (letfn [(metrics-by-id [metrics-fn]
            (reduce (fn [coll x]
                      (assoc coll (:id x) x))
                    {}
                    (metrics-fn vcf-file ref-file :metrics (when metrics (map :id metrics))
                                :use-subsample? use-subsample?)))
          (present-metrics-by-id [base metrics-fn]
            (-> (metrics-by-id metrics-fn)
                (select-keys (keys base))))]
    (let [base-metrics (metrics-by-id im/get-raw-metrics)]
      (->> (merge-with merge
                       base-metrics
                       (present-metrics-by-id base-metrics gemini/get-raw-metrics))
           vals
           (sort-by :id)))))

;; ## API functions

(defn available-metrics
  [file-id & {:keys [rclient]}]
  (let [vcf-file (when file-id (remote/get-file file-id rclient))]
    (concat (im/available-metrics vcf-file)
            (gemini/available-metrics vcf-file))))

(defn plot-ready-metrics
  "Provide metrics for a VCF file ready for plotting and visualization."
  [in-vcf-file & {:keys [metrics rclient]}]
  (let [vcf-file (remote/get-file in-vcf-file rclient)
        ref-file (-> @web-config :ref first :genome)
        plot-metrics (or metrics (available-metrics vcf-file))
        raw-metrics (clean-raw-metrics
                     (combined-raw-metrics vcf-file ref-file plot-metrics false)
                     (map :id plot-metrics))]
    {:filename in-vcf-file
     :created-on (java.util.Date.)
     :metrics (map #(merge % (prepare-plot-metrics (:id %) (get raw-metrics (:id %))))
                   (remove #(nil? (get raw-metrics (:id %))) plot-metrics))}))

(defn get-raw-metrics
  "Retrieve raw metrics values from input VCF."
  [variant-id & {:keys [metrics rclient use-subsample?]}]
  (let [vcf-file (remote/get-file variant-id rclient)
        ref-file (-> @web-config :ref first :genome)]
    (combined-raw-metrics vcf-file ref-file metrics use-subsample?)))
