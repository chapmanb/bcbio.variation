(ns bcbio.variation.api.metrics
  "Provide high level API for accessing variant associated metrics."
  (:import [org.jfree.data.statistics HistogramDataset HistogramType])
  (:use [ordered.map :only [ordered-map]]
        [bcbio.variation.api.file :only [retrieve-file]]
        [bcbio.variation.filter.classify :only [get-vc-attrs]]
        [bcbio.variation.variantcontext :only [get-vcf-header get-vcf-source parse-vcf]]))

(def ^{:doc "Metrics to expose, ranked in order of priority with default min/max values."
       :private true}
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

(defn- get-histogram-bins
  [items n bin-min bin-max]
  "Retrieve values binned into a histogram using JFree Chart."
  (let [ds (doto (HistogramDataset.)
             (.setType HistogramType/RELATIVE_FREQUENCY)
             (.addSeries 0 (double-array items) n bin-min bin-max))]
    {:x (map #(.getXValue ds 0 %) (range (.getItemCount ds 0)))
     :y (map #(.getYValue ds 0 %) (range (.getItemCount ds 0)))}))

(defn- get-raw-metrics
  "Retrieve raw metrics values from input VCF for provided keys."
  [ks vcf-file ref-file]
  (with-open [vcf-source (get-vcf-source vcf-file ref-file)]
    (reduce (fn [coll vc]
              (reduce (fn [inner-coll [key val]]
                        (assoc inner-coll key
                               (cons val (get inner-coll key))))
                      coll (get-vc-attrs vc (keys coll))))
            (zipmap ks (repeat [])) (parse-vcf vcf-source))))

(defn- clean-raw-metrics
  "Remove nil values and empty input metrics."
  [raw]
  (reduce (fn [coll [k vs]]
            (let [clean-vs (remove nil? vs)]
              (if (empty? clean-vs)
                coll
                (assoc coll k clean-vs))))
          {} raw))

(defn- prepare-plot-metrics
  "Bin metrics in preparation for histogram display using predefined min-max boundaries."
  [metric raw]
  (let [bins 20
        [bin-min bin-max] (get expose-metrics metric)
        data (get-histogram-bins raw bins bin-min bin-max)]
    {:vals (:y data)
     :bin-width (- (second (:x data)) (first (:x data)))
     :x-scale {:type :linear
               :domain [bin-min bin-max]}
     :y-scale {:type :linear}}))

(defn plot-ready-metrics
  "Provide metrics for a VCF file ready for plotting and visualization."
  [in-vcf-file ref-file & {:keys [metrics creds cache-dir]}]
  (let [vcf-file (retrieve-file in-vcf-file creds cache-dir)
        plot-metrics (if (nil? metrics) (available-metrics vcf-file) metrics)
        raw-metrics (clean-raw-metrics
                     (get-raw-metrics (map :id plot-metrics) vcf-file ref-file))]
    {:filename vcf-file
     :created-on (java.util.Date.)
     :metrics (map #(merge % (prepare-plot-metrics (:id %) (get raw-metrics (:id %))))
                   (remove #(nil? (get raw-metrics (:id %))) plot-metrics))}))
