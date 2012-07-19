(ns bcbio.variation.filter.trusted
  "Retrieve trusted variants from comparisons based on configured thresholds.
  Allows specification of cases where we should trust variants to pass, such
  as: found in more than two sequencing technologies, or called in 3 aligners,
  or called in 7 out of 8 inputs."
  (:use [bcbio.variation.multiple :only [multiple-overlap-analysis remove-mod-name
                                         prep-cmp-name-lookup]])
  (:require [clojure.string :as string]))

(defn- pairwise-only?
  "Check if a comparison set is only pairwise and not multiple."
  [cmp-names]
  (= 1 (count (set (map (fn [xs] (vec (map remove-mod-name xs))) cmp-names)))))

(defn get-support-vcfs
  "Retrieve supporting VCFs for a set of comparisons and specified support."
  [cmps support config]
  (let [cmps-by-name (if (map? cmps) cmps (prep-cmp-name-lookup cmps))
        support (if (and (not (coll? support)) (pairwise-only? (keys cmps-by-name)))
                  (first (keys cmps-by-name))
                  support)]
    (if (coll? support)
      (zipmap [:true-positives :false-positives]
              (take 2 (-> cmps-by-name (get support) :c-files vals)))
      (let [x (multiple-overlap-analysis cmps-by-name config support)]
        (into {} (map (juxt identity x)
                      [:true-positives :false-positives :target-overlaps]))))))

(defn variant-set-metadata
  "Retrieve metadata associated with overlapping variants from combined set attribute."
  [vc calls]
  (when-let [set-val (get-in vc [:attributes "set"])]
    (let [set-calls (if (= set-val "Intersection")
                      (set (map :name calls))
                      (->> (string/split set-val #"-")
                           (remove #(.startsWith % "filter"))
                           (map #(string/split % #"AND"))
                           (apply concat)
                           set))]
      (reduce (fn [coll x]
                (let [cur-name (string/replace (:name x) "-" "_")]
                  (if-not (contains? set-calls cur-name)
                    coll
                    (reduce (fn [inner [k v]]
                              (assoc inner k (conj (get inner k #{}) v)))
                            coll (assoc (get x :metadata {}) :total cur-name)))))
              {} calls))))

(defn get-trusted-variants
  "Retrieve VCF file of trusted variants based on specific parameters."
  [cmps support params config]
  (when-let [base-vcf (:target-overlaps (get-support-vcfs cmps support config))]
    (println base-vcf)))