(ns bcbio.variation.api.metrics
  "Provide high level API for accessing variant associated metrics."
  (:import [org.jfree.data.statistics HistogramDataset HistogramType])
  (:use [ordered.set :only [ordered-set]]
        [bcbio.variation.filter.classify :only [get-vc-attrs]]
        [bcbio.variation.variantcontext :only [get-vcf-header get-vcf-source parse-vcf]]))

(def ^{:doc "Metrics to expose, ranked in order of priority" :private true}
  expose-metrics
  (ordered-set "QUAL" "DP" "MQ" "QD" "HaplotypeScore"))

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
                                {} (map-indexed vector expose-metrics))]
      (->> (get-vcf-header vcf-file)
           .getMetaData
           (filter #(= "INFO" (.getKey %)))
           (filter #(contains? expose-metrics (.getID %)))
           (map convert-header)
           (concat default-metrics)
           (sort-by #(get metrics-order (:id %)))))))

(defn- get-histogram-bins
  [items n]
  "Retrieve values binned into a histogram using JFree Chart."
  (let [ds (doto (HistogramDataset.)
             (.setType HistogramType/RELATIVE_FREQUENCY)
             (.addSeries 0 (double-array items) n))]
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

(defn- prepare-plot-metrics
  [raw]
  (let [bins 20
        data (get-histogram-bins raw bins)]
    {:vals (:y data)
     :bin-width (- (second (:x data)) (first (:x data)))
     :x-scale {:type :linear
               :domain [(apply min (:x data)) (apply max (:x data))]}
     :y-scale {:type :linear}}))

(defn plot-ready-metrics
  "Provide metrics for a VCF file ready for plotting and visualization."
  [vcf-file ref-file & {:keys [metrics]}]
  (let [plot-metrics (if (nil? metrics) (available-metrics vcf-file) metrics)
        raw-metrics (get-raw-metrics (map :id plot-metrics) vcf-file ref-file)]
    {:filename vcf-file
     :created-on (java.util.Date.)
     :metrics (map #(merge % (prepare-plot-metrics (get raw-metrics (:id %))))
                   plot-metrics)}))
