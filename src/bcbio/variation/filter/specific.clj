(ns bcbio.variation.filter.specific
  "Identify technology or caller specific variants from multiple combined callsets."
  (:use [bcbio.variation.filter.trusted :only [variant-set-metadata]]))

(defn- get-specific
  [data kw-want kw-cmp]
  (when (and (= 1 (count (kw-want data)))
             (> (count (kw-cmp data)) 1))
    (first (kw-want data))))

(defn get-x-specific-designation
  "Check if a variant is specific to a caller or method."
  [vc calls]
  (let [data (variant-set-metadata vc calls)]
    (reduce (fn [coll [kw-want kw-cmp]]
              (if-let [x (get-specific data kw-want kw-cmp)]
                (assoc coll kw-want x)
                coll))
            {} [[:caller :technology] [:technology :caller]])))