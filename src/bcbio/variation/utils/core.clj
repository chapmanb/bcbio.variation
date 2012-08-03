(ns bcbio.variation.utils.core
  (:require [bcbio.variation.utils.popfreq :as popfreq]
            [bcbio.variation.utils.callsummary :as callsummary]
            [bcbio.variation.utils.summarize :as summarize]))

(defn -main [cur-type & args]
  (apply (case (keyword cur-type)
           :popfreq popfreq/annotate-with-popfreq
           :callsummary callsummary/annotate-with-callsummary
           :summarize summarize/vcf-to-table-config)
         args))