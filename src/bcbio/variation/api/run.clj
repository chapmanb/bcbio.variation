(ns bcbio.variation.api.run
  "High level API to run analyses."
  (:use [bcbio.variation.filter :only [variant-filter jexl-filters-from-map]]
        [bcbio.variation.api.shared :only [web-config]])
  (:require [fs.core :as fs]
            [bcbio.variation.api.file :as file-api]))

(defmulti do-analysis
  "Run analysis on provided inputs, dispatching on analysis type"
  (fn [atype params creds] (keyword atype)))

(defmethod do-analysis :filter
  ^{:doc "Filter an input file according to specified metrics.
          params:
            - filename: The file to process
            - metrics: A map of filters, with metrics names as keys
              and [min max] as values. Long term we could dispatch on
              different value types for categorical data."}
  [atype params creds]
  (let [genome (-> @web-config :ref first)
        in-file (file-api/retrieve-file (:filename params) creds
                                        (get-in @web-config [:dir :cache]))
        filter-file (variant-filter in-file
                                    (jexl-filters-from-map (:metrics params))
                                    (:genome genome))
        local-out-dir (fs/file (fs/parent in-file) (name atype))
        remote-dir (file-api/put-files [filter-file] (:filename params) (name atype)
                                       creds)]
    (when-not (fs/exists? local-out-dir)
      (fs/mkdirs local-out-dir))
    (doseq [ext ["" ".idx"]]
      (fs/rename (str filter-file ext) (str (fs/file local-out-dir (fs/base-name filter-file)) ext)))
    (file-api/get-files :vcf creds [remote-dir])))

(defmethod do-analysis :score
  ^{:doc "Run comparison and scoring analysis on provided input files.
          params:
            - gs-variant-file: Input variant file, from GenomeSpace
            - gs-region-file: Optional BED file of regions to score on. From GenomeSpace.
            - comparison-genome: Name of genome to compare against. Used
              to look up comparison details in configuration file."}
  [_ params creds]
  (println params))