(ns bcbio.variation.utils.core
  (:require [bcbio.variation.normalize :as normalize]
            [bcbio.variation.utils.callsummary :as callsummary]
            ;[bcbio.variation.utils.gms :as gms]
            [bcbio.variation.utils.illumina :as illumina]
            [bcbio.variation.utils.popfreq :as popfreq]
            [bcbio.variation.utils.quickcompare :as qcmp]
            [bcbio.variation.utils.summarize :as summarize]
            [bcbio.variation.utils.svmerge :as svmerge]))

(def ^{:private true} progs
  {:callsummary callsummary/annotate-with-callsummary
                                        ;:gms gms/prepare-gms-vcfs-from-config
   :illumina illumina/cl-entry
   :popfreq popfreq/annotate-with-popfreq
   :sort-vcf normalize/prep-vcf
   :summarize summarize/vcf-to-table-config
   :svmerge svmerge/into-calls
   :quickcompare qcmp/-main})

(defn -main [& raw-args]
  (if (or (empty? raw-args) (not (contains? progs (keyword (first raw-args)))))
    (do
      (println "variant-utils requires at least one argument for the sub-command to run:")
      (doseq [prog (keys progs)]
        (println (name prog))))
    (apply (get progs (keyword (first raw-args)))
           (rest raw-args))))
