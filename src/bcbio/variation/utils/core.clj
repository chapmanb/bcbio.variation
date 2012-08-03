(ns bcbio.variation.utils.core
  (:require [bcbio.variation.utils.callsummary :as callsummary]
            [bcbio.variation.utils.gms :as gms]
            [bcbio.variation.utils.popfreq :as popfreq]
            [bcbio.variation.utils.summarize :as summarize]))

(defn -main [cur-type & args]
  (apply (case (keyword cur-type)
           :callsummary callsummary/annotate-with-callsummary
           :gms gms/prepare-gms-vcfs-from-config
           :popfreq popfreq/annotate-with-popfreq
           :summarize summarize/vcf-to-table-config)
         args))