(ns bcbio.variation.api.run
  "High level API to run analyses."
  (:use [bcbio.variation.filter :only [variant-filter jexl-filters-from-map]]
        [bcbio.variation.api.shared :only [web-config]])
  (:require [clojure.string :as string]
            [fs.core :as fs]
            [bcbio.variation.api.file :as file-api]
            [bcbio.variation.remote.core :as remote]
            [bcbio.variation.workflow.xprize :as xprize]))

(defmulti do-analysis
  "Run analysis on provided inputs, dispatching on analysis type"
  (fn [atype params rclient] (keyword atype)))

(defn- run-filter
  "Run filtering, pushing results to remote file store.
   Returns list of output files following pushing the filter."
  [atype params rclient]
  (let [ref-file (-> @web-config :ref first :genome)
        in-file (remote/get-file (:filename params) rclient)
        filter-file (variant-filter in-file
                                    (jexl-filters-from-map (:metrics params))
                                    ref-file)
        local-out-dir (fs/file (fs/parent in-file) (name atype))
        remote-dir (remote/put-file rclient filter-file {:input-file (:filename params)
                                                         :tag (name atype)
                                                         :file-type :vcf})]
    (when-not (fs/exists? local-out-dir)
      (fs/mkdirs local-out-dir))
    (doseq [ext ["" ".idx"]]
      (fs/rename (str filter-file ext) (str (fs/file local-out-dir (fs/base-name filter-file)) ext)))
    (remote/list-files rclient remote-dir :vcf)))

(defmethod do-analysis :filter
  ^{:doc "Filter an input file according to specified metrics.
          params:
            - filename: The file to process
            - metrics: A map of filters, with metrics names as keys
              and [min max] as values. Long term we could dispatch on
              different value types for categorical data."}
  [atype params rclient]
  {:runner (future (run-filter atype params rclient))})

(defn- prep-xprize-files
  "Prepare X prize input files, potentially pulling from remote directories."
  [work-info rclient]
  (letfn [(get-remote-files [work-info]
            (reduce (fn [coll kw]
                      (assoc coll kw (remote/get-file (get coll kw) rclient
                                                      :out-dir (:dir work-info))))
                    work-info [:variant-file :region-file]))]
    (-> work-info
        (assoc :orig-variant-file (:variant-file work-info))
        get-remote-files)))

(defn- upload-xprize-files
  "Upload X Prize results files back to remote directories."
  [{:keys [work-info comparison]} rclient host-info]
  (when-not (nil? (:conn rclient))
    (doseq [x (map #(get-in comparison [:c-files %])
                   [:concordant :discordant :discordant-missing :phasing-error :summary])]
      (remote/put-file rclient x {:dbkey "hg_g1k_v37"
                                  :file-type "hg19"
                                  :input-file (:orig-variant-file work-info)
                                  :host-info host-info
                                  :tag "xprize"})))
   comparison)

(defmethod do-analysis :xprize
  ^{:doc "Run X Prize comparison and scoring analysis on provided input files.
          params:
            - variant-file: Input variant file, in VCF format, to compare.
            - region-file: Optional BED file of regions to score on.
            - comparison-genome: Name of genome to compare against. Used
              to look up comparison details in configuration file.
            - host-info: Host information for providing callbacks to a local server."}
  [atype params rclient]
  (let [work-info (xprize/prep-scoring params @web-config)]
    {:runner (future (-> work-info
                         (prep-xprize-files rclient)
                         (xprize/run-scoring-analysis rclient
                                                      (:host-info params) @web-config)
                         (upload-xprize-files rclient (:host-info params))))
     :work-info work-info}))