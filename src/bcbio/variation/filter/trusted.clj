(ns bcbio.variation.filter.trusted
  "Retrieve trusted variants from comparisons based on configured thresholds.
  Allows specification of cases where we should trust variants to pass, such
  as: found in more than two sequencing technologies, or called in 3 aligners,
  or called in 7 out of 8 inputs."
  (:use [bcbio.variation.multiple :only [multiple-overlap-analysis remove-mod-name]]))

(defn- pairwise-only?
  "Check if a comparison set is only pairwise and not multiple."
  [cmp-names]
  (= 1 (count (set (map (fn [xs] (vec (map remove-mod-name xs))) cmp-names)))))

(defn get-support-vcfs
  "Retrieve supporting VCFs for a set of comparisons and specified support."
  [cmps-by-name support config]
  (let [support (if (and (not (coll? support)) (pairwise-only? (keys cmps-by-name)))
                  (first (keys cmps-by-name))
                  support)]
    (if (coll? support)
      (zipmap [:true-positives :false-positives]
              (take 2 (-> cmps-by-name (get support) :c-files vals)))
      (let [x (multiple-overlap-analysis cmps-by-name config support)]
        (into {} (map (comp identity x)
                      [:true-positives :false-positives :target-overlaps]))))))

(defn get-trusted-variants
  "Retrieve VCF file of trusted variants based on specific parameters."
  [cmps support params config]
  (when-let [base-vcf (:target-overlaps (get-support-vcfs cmps support config))]
    (println base-vcf)))