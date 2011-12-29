;; Parse and provide detailed information from GATKReport output

(ns bcbio.variation.report
  (:import [org.broadinstitute.sting.gatk.report GATKReport]))

(defn concordance-report-metrics [sample in-file]
  (letfn [(sample-in-row? [x]
            (and (= (:row x) sample)
                 (= (:Sample x) sample)
                 (= (:Novelty x) "all")))]
    (let [table (-> (GATKReport. in-file) (.getTable "GenotypeConcordance.simplifiedStats"))
          cols (rest (.getColumns table))
          headers (map #(-> % (.getColumnName) keyword) cols)]
      (filter sample-in-row?
              (for [i (range (count (.values (first cols))))]
                (zipmap headers
                        (map #(nth (vec (.values %)) i) cols)))))))

(defn write-concordance-metrics [metrics wrtr]
  (.write wrtr (format "Overall genotype concordance   %s\n"
                    (:percent_overall_genotype_concordance metrics)))
  (.write wrtr (format "Non-reference discrepancy rate %s\n"
                       (:percent_non_reference_discrepancy_rate metrics)))
  (.write wrtr (format "Non-reference sensitivity      %s\n"
                       (:percent_non_reference_sensitivity metrics))))
