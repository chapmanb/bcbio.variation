(ns bcbio.variation.utils.core
  (:require [bcbio.variation.utils.callsummary :as callsummary]
            ;[bcbio.variation.utils.gms :as gms]
            [bcbio.variation.utils.illumina :as illumina]
            [bcbio.variation.utils.popfreq :as popfreq]
            [bcbio.variation.utils.quickcompare :as qcmp]
            [bcbio.variation.utils.summarize :as summarize]
            [bcbio.variation.utils.svmerge :as svmerge]))

(defn -main [cur-type & args]
  (apply (case (keyword cur-type)
           :callsummary callsummary/annotate-with-callsummary
           ;:gms gms/prepare-gms-vcfs-from-config
           :illumina illumina/cl-entry
           :popfreq popfreq/annotate-with-popfreq
           :summarize summarize/vcf-to-table-config
           :svmerge svmerge/into-calls
           :quickcompare qcmp/-main)
         args))
