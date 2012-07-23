(ns bcbio.variation.filter.trusted
  "Retrieve trusted variants from comparisons based on configured thresholds.
  Allows specification of cases where we should trust variants to pass, such
  as: found in more than two sequencing technologies, or called in 3 aligners,
  or called in 7 out of 8 inputs."
  (:use [bcbio.variation.multiple :only [multiple-overlap-analysis remove-mod-name
                                         prep-cmp-name-lookup]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template
                                               get-vcf-source]])
  (:require [clojure.string :as string]
            [bcbio.run.itx :as itx]))

(defn- pairwise-only?
  "Check if a comparison set is only pairwise and not multiple."
  [cmp-names]
  (= 1 (count (set (map (fn [xs] (vec (map remove-mod-name xs))) cmp-names)))))

(defn get-support-vcfs
  "Retrieve supporting VCFs for a set of comparisons and specified support."
  [cmps support config & {:keys [remove-mods?]}]
  (let [cmps-by-name (if (map? cmps) cmps (prep-cmp-name-lookup cmps :remove-mods? remove-mods?))
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

(defn is-trusted-variant?
  "Determine if we trust a variant based on specified trust parameters.
  The params specify required counts for inclusion. For instance:
  {:total 4 :technology 3 :caller 2} includes variants located in 4 total calls
  or in three different technologies or in 2 different callers.
  It can also handle percentages for required inputs:
  {:total 1.0 :technology 0.75}"
  [vc params calls]
  (letfn [(collapse-md-by-type [calls]
            (reduce (fn [coll [k v]]
                      (assoc coll k (conj (get coll k #{}) v)))
                    {:total (set (map :name calls))}
                    (mapcat :metadata calls)))
          (calc-md-counts [calls]
            (reduce (fn [coll [k v]]
                      (assoc coll k (count v)))
                    {}
                    (collapse-md-by-type calls)))
          (param-passes? [metadata md-counts [k v]]
            (let [n (count (get metadata k []))]
              (if (> v 1)
                (>= n v)
                (>= (/ n (get md-counts k)) v))))]
    (some (partial param-passes? (variant-set-metadata vc calls) (calc-md-counts calls))
          params)))

(defn get-trusted-variants
  "Retrieve VCF file of trusted variants based on specific parameters."
  [cmps support params exp config]
  (when-let [base-vcf (:target-overlaps
                       (get-support-vcfs cmps (if (coll? support) (first support) support)
                                         config :remove-mods? true))]
    (let [out-file (itx/add-file-part base-vcf "trusted")]
      (when (itx/needs-run? out-file)
        (with-open [base-vcf-s (get-vcf-source base-vcf (:ref exp))]
          (write-vcf-w-template base-vcf {:out out-file}
                                (->> (parse-vcf base-vcf-s)
                                     (filter #(is-trusted-variant? % params (:calls exp)))
                                     (map :vc))
                                (:ref exp))))
      out-file)))