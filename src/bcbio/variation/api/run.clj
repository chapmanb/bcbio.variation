(ns bcbio.variation.api.run
  "High level API to run analyses."
  (:use [bcbio.variation.filter :only [category-variant-filter]]
        [bcbio.variation.api.shared :only [web-config]])
  (:require [clojure.string :as string]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.variation.api.file :as file-api]
            [bcbio.variation.remote.core :as remote]
            [bcbio.variation.workflow.xprize :as xprize]))

(defmulti do-analysis
  "Run analysis on provided inputs, dispatching on analysis type"
  (fn [atype params rclient] (keyword atype)))

(defn- to-prep-file [x]
  (itx/add-file-part x "prep"))

(defn- from-prep-file [x]
  (itx/remove-file-part x "prep"))

(defn- run-filter
  "Run filtering, pushing results to remote file store.
   Returns list of output files following pushing the filter."
  [atype params rclient]
  (let [ref-file (-> @web-config :ref first :genome)
        in-file (to-prep-file (remote/get-file (:filename params) rclient))
        filter-file (category-variant-filter in-file (:metrics params) ref-file :remove? true)
        local-out-dir (fs/file (fs/parent in-file) (name atype))
        final-filter-file (str (fs/file local-out-dir (from-prep-file (fs/base-name filter-file))))]
    (when-not (fs/exists? local-out-dir)
      (fs/mkdirs local-out-dir))
    (doseq [ext ["" ".idx"]]
      (fs/rename (str filter-file ext) (str final-filter-file ext)))
    (doseq [ext ["-prep-gemini.db" "-prep-metrics.db" "-prep.vcf"]]
      (let [idx-name (str (fs/file local-out-dir (fs/name final-filter-file)) ext)]
        (when (fs/exists? idx-name)
          (fs/delete idx-name))))
    (remote/put-file rclient final-filter-file
                     {:input-file (:filename params)
                      :expose-fn (:expose-fn params)
                      :tag (name atype)
                      :file-type :vcf})
    ;; Do not need to re-prep this file since it's a subset of a prep
    (fs/touch final-filter-file)
    (fs/copy final-filter-file (to-prep-file final-filter-file))
    (future (file-api/pre-fetch-remotes rclient))
    final-filter-file))

(defmethod do-analysis :filter
  ^{:doc "Filter an input file according to specified metrics.
          params:
            - filename: The file to process
            - metrics: A map of filters, with metrics names as keys
              and either ranges ([min max]) or categories as values."}
  [atype params rclient]
  {:runner (future (to-prep-file
                    (if (empty? (:metrics params))
                      (remote/get-file (:filename params) rclient)
                      (run-filter atype params rclient))))})

(defn- prep-xprize-files
  "Prepare X prize input files, potentially pulling from remote directories."
  [work-info rclient]
  (letfn [(get-remote-files [work-info]
            (reduce (fn [coll kw]
                      (assoc coll kw (when-let [f (get coll kw)]
                                       (remote/get-file f rclient :out-dir (:dir work-info)))))
                    work-info [:variant-file :region-file]))]
    (-> work-info
        (assoc :orig-variant-file (:variant-file work-info))
        get-remote-files)))

(defn- upload-xprize-files
  "Upload X Prize results files back to remote directories."
  [{:keys [work-info comparison]} rclient params]
  (when-not (nil? (:conn rclient))
    (when-let [remote-input (:orig-variant-file work-info)]
      (doseq [x (map #(get-in comparison [:c-files %])
                     [:concordant :discordant :discordant-missing :phasing-error :summary])]
        (let [ftype (cond
                     (.endsWith x ".vcf") :vcf
                     :else :tabular)]
          (remote/put-file rclient x {:dbkey :hg19
                                      :file-type ftype
                                      :input-file remote-input
                                      :expose-fn (:expose-fn params)
                                      :tag "xprize"})))))
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
                         (xprize/run-scoring-analysis rclient @web-config)
                         (upload-xprize-files rclient params)))
     :work-info work-info}))