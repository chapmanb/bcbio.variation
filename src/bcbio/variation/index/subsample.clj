(ns bcbio.variation.index.subsample
  "Provide rapid subsampling capabilities for indexed retrieval of metrics.
  Used on index preparation to provide a representative subset of large
  datasets based on assessment metrics. Sub-sampled metrics allow
  interactive visualizations of large data."
  (:require [clj-ml.data :as mldata]
            [clj-ml.clusterers :as mlclust]))

(defn- final-subsample-ids
  "Return ids of subsampled metrics with a single representative from each cluster."
  [xs clusters]
  (let [final-is (-> (reduce (fn [coll [i cluster]]
                               {:seen (conj (:seen coll) cluster)
                                :want (if (contains? (:seen coll) cluster)
                                        (:want coll)
                                        (conj (:want coll) i))})
                             {:seen #{} :want []} (map-indexed vector clusters))
                     :want
                     set)]
    (remove nil?
            (map-indexed (fn [i x] (when (contains? final-is i) x)) xs))))

(defn subsample-by-cluster
  "Subsample a set of metrics using clustering. Returns ids of
  representative items from each cluster."
  [metrics params]
  (letfn [(get-attrs [attrs x] (map #(get x %) attrs))]
    (let [clusterer (mlclust/make-clusterer (keyword (get-in params [:subsample :method]))
                                            {:number-clusters (get-in params [:subsample :count])})
          attrs (-> metrics first (dissoc :id) keys)
          ds (mldata/make-dataset "ds" attrs
                                  (map (partial get-attrs attrs) metrics))]
      (mlclust/clusterer-build clusterer ds)
      (->> (mlclust/clusterer-cluster clusterer ds)
           mldata/dataset-seq
           (map mldata/instance-get-class)
           (final-subsample-ids (map :id metrics))))))