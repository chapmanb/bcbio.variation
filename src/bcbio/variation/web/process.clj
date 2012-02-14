(ns bcbio.variation.web.process
  "Run scoring analysis, handling preparation of input files and run configuration."
  (:import [java.util.UUID])
  (:use [clojure.java.io]
        [ring.util.response :only (redirect)]
        [bcbio.variation.compare :only (variant-comparison-from-config)])
  (:require [fs.core :as fs]
            [clj-yaml.core :as yaml]))

(defn create-work-config
  "Create configuration for processing inputs using references supplied in config."
  [in-files work-info config]
  (if-not (fs/exists? (:dir work-info))
    (fs/mkdirs (:dir work-info)))
  (let [config-file (str (fs/file (:dir work-info) "process.yaml"))]
    (->> {:outdir (str (fs/file (:dir work-info) "grading"))
          :outdir-prep (str (fs/file (:dir work-info) "grading" "prep"))
          :experiments [{:sample (-> config :ref :sample)
                         :ref (-> config :ref :genome)
                         :intervals (-> config :ref :intervals)
                         :calls [{:name "reference"
                                  :file (-> config :ref :variants)}
                                 {:name "contestant"
                                  :file (:variant-file in-files)
                                  :intervals (:region-file in-files)}]}]}
         yaml/generate-string
         (spit config-file))
    config-file))

(defn prep-and-run-scoring
  "Prepare output directory and run scoring analysis using form inputs."
  [request]
  (letfn [(prep-tmp-dir [request]
            (let [tmp-dir (-> request :config :dir :work)
                  work-id (str (java.util.UUID/randomUUID))
                  cur-dir (fs/file tmp-dir work-id)]
              (fs/mkdirs cur-dir)
              {:id work-id :dir (str cur-dir)}))
          (download-file [tmp-dir request name]
            (let [cur-param (-> request :params (get name))
                  out-file (fs/file tmp-dir (:filename cur-param))]
              [(keyword name) (if (> (:size cur-param) 0)
                                (do
                                  (copy (:tempfile cur-param) out-file)
                                  (str out-file)))]))]
    (let [work-info (prep-tmp-dir request)
          in-files (into {} (map (partial download-file (:dir work-info) request)
                                 ["variant-file" "region-file"]))
          process-config (create-work-config in-files work-info (:config request))]
      (variant-comparison-from-config process-config)))
  (redirect "/"))
