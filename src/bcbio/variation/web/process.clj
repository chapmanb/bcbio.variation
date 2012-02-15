(ns bcbio.variation.web.process
  "Run scoring analysis, handling preparation of input files and run configuration."
  (:import [java.util.UUID])
  (:use [clojure.java.io]
        [ring.util.response :only (redirect)]
        [bcbio.variation.compare :only (variant-comparison-from-config)]
        [bcbio.variation.report :only (prep-scoring-table)])
  (:require [fs.core :as fs]
            [clj-yaml.core :as yaml]
            [net.cgrand.enlive-html :as html]
            [doric.core :as doric]))

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
                                  :file (if-let [x (:variant-file in-files)]
                                          x
                                          (-> config :ref :default-compare))
                                  :intervals (if-let [x (:region-file in-files)]
                                               x
                                               (-> config :ref :intervals))}]}]}
         yaml/generate-string
         (spit config-file))
    config-file))

(defn html-scoring-summary
  "Generate a summary of scoring results for display."
  [comparisons]
  {:pre [(= 1 (count comparisons))]}
  (let [scoring-table (prep-scoring-table (-> comparisons first :metrics))]
    (-> (str (doric/table ^{:format doric/html} [:metric :value] scoring-table))
     java.io.StringReader.
     html/html-resource
     (html/transform [:table] (html/set-attr :class "table table-condensed")))))

(defn update-main-html
  "Update main page HTML with replaced main panel content."
  [request new-content]
  (apply str (html/emit*
              (html/transform (html/html-resource (fs/file (-> request :config :dir :html-root)
                                                           "index.html"))
                              [:div#main-content]
                              (html/content new-content)))))

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
          process-config (create-work-config in-files work-info (:config request))
          comparisons (variant-comparison-from-config process-config)]
      (update-main-html request
                        (html-scoring-summary comparisons)))))
