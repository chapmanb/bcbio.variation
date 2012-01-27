;; Parse and provide detailed information from GATKReport output

(ns bcbio.variation.report
  (:import [org.broadinstitute.sting.gatk.report GATKReport])
  (:use [ordered.map :only [ordered-map]])
  (:require [doric.core :as doric]))

(defn concordance-report-metrics [sample in-file]
  "Retrieve high level concordance metrics from GATK VariantEval report."
  (letfn [(sample-in-row? [x]
            (and (= (:row x) sample)
                 (= (:Sample x) sample)
                 (= (:Novelty x) "all")
                 (= (:Filter x) "called")))]
    (let [table (-> (GATKReport. in-file) (.getTable "GenotypeConcordance.simplifiedStats"))
          cols (rest (.getColumns table))
          headers (map #(-> % (.getColumnName) keyword) cols)]
      (filter sample-in-row?
              (for [i (range (count (.values (first cols))))]
                (zipmap headers
                        (map #(nth (vec (.values %)) i) cols)))))))

(defn write-concordance-metrics [metrics wrtr]
  (let [to-write {:genotype_concordance "Overall genotype concordance"
                  :nonref_discrepency "Non-reference discrepancy rate"
                  :nonref_sensitivity "Non-reference sensitivity"
                  :concordant "Concordant count"
                  :nonref_concordant "Non-reference concordant count"
                  :discordant1 (str "Discordant: " (:call1 metrics))
                  :discordant2 (str "Discordant: " (:call2 metrics))}]
    (.write wrtr (str (doric/table [:metric :value]
                                   (map (fn [[k v]] {:metric v :value (get metrics k)})
                                        to-write))
                      "\n"))))
