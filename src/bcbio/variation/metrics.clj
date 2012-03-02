(ns bcbio.variation.metrics
   "Accumulate and analyze metrics associated with each variant.
   This provides summaries intended to identify characteristic
   metrics to use for filtering."
  (:use [clojure.java.io]
        [bcbio.variation.variantcontext :only [parse-vcf get-vcf-source]])
  (:require [incanter.stats :as istats]
            [doric.core :as doric]))

;; ## Convenience functions

(defn passes-filter? [vc]
  (= (count (:filters vc)) 0))

;; ## Summary metrics
;; Provide a summary-style presentation of distribution of metrics values.

(defn- to-float [x]
  (try
    (Float/parseFloat x)
    (catch Exception e x)))

(def header [{:name :metric}
             {:name :count}
             {:name :min :format #(format "%.2f" %)}
             {:name :pct25 :format #(format "%.2f" %)}
             {:name :median :format #(format "%.2f" %)}
             {:name :pct75 :format #(format "%.2f" %)}
             {:name :max :format #(format "%.2f" %)}])

(defn summary-stats [key vals]
  "Provide summary statistics on a list of values."
  (zipmap (map :name header)
          (concat [key (count vals)]
                  (istats/quantile vals))))

(defn- raw-vcf-stats [vcf-file]
  "Accumulate raw statistics associated with variant calls from input VCF."
  (letfn [(collect-attributes [collect [k v]]
            (if (number? (to-float v))
              (assoc collect k (cons (to-float v) (get collect k [])))
              collect))
          (collect-vc [collect vc]
            (assoc (reduce collect-attributes collect (:attributes vc))
              "QUAL" (cons (-> vc :genotypes first :qual)
                           (get collect "QUAL" []))))]
    (with-open [vcf-source (get-vcf-source vcf-file)]
      (reduce collect-vc {} (filter passes-filter? (parse-vcf vcf-source))))))

(defn vcf-stats [vcf-file]
  "Collect summary statistics associated with variant calls."
  (let [raw-stats (raw-vcf-stats vcf-file)]
    (map #(apply summary-stats %) (sort-by first raw-stats))))

(defn write-summary-table [stats & {:keys [wrtr]
                                    :or {wrtr (writer System/out)}}]
  (.write wrtr (str (doric/table header stats) "\n")))

;; ## Classify
;; Provide metrics for files in preparation for automated
;; classification.

(defn get-vcf-classifier-metrics
  "Collect metrics into tables ready to feed into classification algorithms."
  [vcf-file]
  (with-open [vcf-source (get-vcf-source vcf-file)]))
