(ns bcbio.variation.utils.core
  (:require [clojure.tools.cli :refer [cli]]
            [bcbio.variation.normalize :as normalize]
            [bcbio.variation.utils.callsummary :as callsummary]
            [bcbio.variation.utils.comparetwo :as comparetwo]
            ;[bcbio.variation.utils.gms :as gms]
            [bcbio.variation.utils.illumina :as illumina]
            [bcbio.variation.utils.popfreq :as popfreq]
            [bcbio.variation.utils.quickcompare :as qcmp]
            [bcbio.variation.utils.summarize :as summarize]
            [bcbio.variation.utils.svmerge :as svmerge]))

(defn- sort-vcf
  "Command line interface to providing preparation and sorting."
  [& args]
  (let [[options [vcf-file ref-file] _]
        (cli args
             ["-s" "--sortpos" "Sort by position" :flag true])]
    (normalize/prep-vcf vcf-file ref-file nil
                        :config {:prep-sort-pos (:sortpos options)
                                 :prep-sv-genotype true
                                 :remove-refcalls false
                                 :prep-org :default})))

(def ^{:private true} progs
  {:callsummary callsummary/annotate-with-callsummary
                                        ;:gms gms/prepare-gms-vcfs-from-config
   :comparetwo comparetwo/cl-entry
   :illumina illumina/cl-entry
   :popfreq popfreq/annotate-with-popfreq
   :sort-vcf sort-vcf
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
