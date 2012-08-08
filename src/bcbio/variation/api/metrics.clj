(ns bcbio.variation.api.metrics
  "Provide high level API for accessing variant associated metrics."
  (:import [org.jfree.data.statistics HistogramDataset HistogramType])
  (:use [bcbio.variation.api.file :only [retrieve-file]]
        [bcbio.variation.filter.classify :only [get-vc-attrs]]
        [bcbio.variation.variantcontext :only [get-vcf-iterator parse-vcf]])
  (:require [bcbio.variation.index.metrics :as im]))

;; ## Helper functions

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

(defn- prepare-plot-metrics
  "Bin metrics in preparation for histogram display using predefined min-max boundaries."
  [metric raw]
  (let [bins 20
        [bin-min bin-max] (:range (get im/expose-metrics metric))
        data (get-histogram-bins raw bins bin-min bin-max)]
    {:vals (:y data)
     :bin-width (- (second (:x data)) (first (:x data)))
     :x-scale {:type :linear
               :domain [bin-min bin-max]}
     :y-scale {:type :linear}}))

;; ## Available API functions

(defn available-metrics
  [file-id & {:keys [creds cache-dir]}]
  (let [vcf-file (retrieve-file file-id creds cache-dir)]
    (im/available-metrics vcf-file)))

(defn plot-ready-metrics
  "Provide metrics for a VCF file ready for plotting and visualization."
  [in-vcf-file ref-file & {:keys [metrics creds cache-dir]}]
  (let [vcf-file (retrieve-file in-vcf-file creds cache-dir)
        plot-metrics (if (nil? metrics) (im/available-metrics vcf-file) metrics)
        raw-metrics (clean-raw-metrics
                     (im/get-raw-metrics vcf-file ref-file :metrics (map :id plot-metrics))
                     (map :id plot-metrics))]
    {:filename in-vcf-file
     :created-on (java.util.Date.)
     :metrics (map #(merge % (prepare-plot-metrics (:id %) (get raw-metrics (:id %))))
                   (remove #(nil? (get raw-metrics (:id %))) plot-metrics))}))

(defn get-raw-metrics
  "Retrieve raw metrics values from input VCF."
  [variant-id ref-file & {:keys [metrics creds cache-dir]}]
  (let [vcf-file (retrieve-file variant-id creds cache-dir)]
    (im/get-raw-metrics vcf-file ref-file :metrics (when metrics
                                                     (map :id metrics)))))
