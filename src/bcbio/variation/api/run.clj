(ns bcbio.variation.api.run
  "High level API to run analyses."
  (:use [bcbio.variation.filter :only [variant-filter jexl-filters-from-map]]
        [bcbio.variation.api.shared :only [web-config]])
  (:require [clojure.string :as string]
            [fs.core :as fs]
            [bcbio.variation.api.file :as file-api]
            [bcbio.variation.remote.core :as remote]))

(defmulti do-analysis
  "Run analysis on provided inputs, dispatching on analysis type"
  (fn [atype params rclient] (keyword atype)))

(defmethod do-analysis :filter
  ^{:doc "Filter an input file according to specified metrics.
          params:
            - filename: The file to process
            - metrics: A map of filters, with metrics names as keys
              and [min max] as values. Long term we could dispatch on
              different value types for categorical data."}
  [atype params rclient]
  (let [ref-file (-> @web-config :ref first :genome)
        in-file (remote/get-file (:filename params) rclient)
        filter-file (variant-filter in-file
                                    (jexl-filters-from-map (:metrics params))
                                    ref-file)
        local-out-dir (fs/file (fs/parent in-file) (name atype))
        remote-dir (str (fs/file (fs/parent (last (string/split (:filename params) #":" 2)))
                              (name atype)))]
    (remote/put-file rclient remote-dir filter-file {})
    (when-not (fs/exists? local-out-dir)
      (fs/mkdirs local-out-dir))
    (doseq [ext ["" ".idx"]]
      (fs/rename (str filter-file ext) (str (fs/file local-out-dir (fs/base-name filter-file)) ext)))
    (remote/list-files rclient remote-dir :vcf)))

(defmethod do-analysis :score
  ^{:doc "Run comparison and scoring analysis on provided input files.
          params:
            - gs-variant-file: Input variant file, from GenomeSpace
            - gs-region-file: Optional BED file of regions to score on. From GenomeSpace.
            - comparison-genome: Name of genome to compare against. Used
              to look up comparison details in configuration file."}
  [_ params rclient]
  (println params))