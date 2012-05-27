(ns bcbio.variation.web.process
  "Run scoring analysis, handling preparation of input files and run configuration."
  (:import [java.util.UUID])
  (:use [clojure.java.io]
        [bcbio.variation.compare :only [variant-comparison-from-config]]
        [bcbio.variation.report :only [prep-scoring-table]]
        [bcbio.variation.web.shared :only [web-config]])
  (:require [fs.core :as fs]
            [clj-yaml.core :as yaml]
            [noir.session :as session]
            [noir.response :as response]
            [net.cgrand.enlive-html :as html]
            [doric.core :as doric]))

;; ## Run scoring based on inputs from web or API

(defn create-work-config
  "Create configuration for processing inputs using references supplied in config."
  [in-files work-info config]
  (if-not (fs/exists? (:dir work-info))
    (fs/mkdirs (:dir work-info)))
  (let [config-file (str (fs/file (:dir work-info) "process.yaml"))]
    (->> {:dir {:out (str (fs/file (:dir work-info) "grading"))
                :prep (str (fs/file (:dir work-info) "grading" "prep"))}
          :experiments [{:sample (-> config :ref :sample)
                         :ref (-> config :ref :genome)
                         :intervals (-> config :ref :intervals)
                         :approach "grade"
                         :calls [{:name "reference"
                                  :file (-> config :ref :variants)
                                  :remove-refcalls true
                                  :refcalls false}
                                 {:name "contestant"
                                  :prep true
                                  :preclean true
                                  :remove-refcalls true
                                  :file (if-let [x (:variant-file in-files)]
                                          x
                                          (-> config :ref :default-compare))
                                  :intervals (if-let [x (:region-file in-files)]
                                               x
                                               (-> config :ref :intervals))}]}]}
         yaml/generate-string
         (spit config-file))
    config-file))

(defn html-summary-table
  "Generate a summary table of scoring results."
  [comparisons]
  {:pre [(= 1 (count comparisons))]}
  (let [scoring-table (prep-scoring-table (-> comparisons first :metrics))]
    (-> (str (doric/table ^{:format doric/html} [:metric :value] scoring-table))
     java.io.StringReader.
     html/html-resource
     (html/transform [:table] (html/set-attr :class "table table-condensed")))))

(defn html-scoring-summary
  "Generate summary of scoring results for display."
  [comparisons]
  (let [template-dir (get-in @web-config [:dir :template])
        sum-table (html-summary-table comparisons)]
    (html/transform (html/html-resource (fs/file template-dir "scoresum.html"))
                    [:div#score-table]
                    (html/content sum-table))))

(defn scoring-html
  "Update main page HTML with content for scoring."
  []
  (let [html-dir (get-in @web-config [:dir :html-root])
        template-dir (get-in @web-config [:dir :template])]
    (apply str (html/emit*
                (html/transform (html/html-resource (fs/file html-dir "index.html"))
                                [:div#main-content]
                                (-> (fs/file template-dir "score.html")
                                    html/html-resource
                                    (html/select [:body])
                                    first
                                    html/content))))))

(defn run-scoring
  "Run scoring analysis from details provided in current session."
  []
  (let [process-config (create-work-config (session/get :in-files)
                                           (session/get :work-info)
                                           @web-config)
        comparisons (variant-comparison-from-config process-config)]
    (apply str (html/emit*
                (html-scoring-summary comparisons)))))

(defn prep-scoring
  "Download form-supplied input files and prep directory for scoring analysis."
  [params]
  (letfn [(prep-tmp-dir []
            (let [tmp-dir (get-in @web-config [:dir :work])
                  work-id (str (java.util.UUID/randomUUID))
                  cur-dir (fs/file tmp-dir work-id)]
              (fs/mkdirs cur-dir)
              {:id work-id :dir (str cur-dir)}))
          (download-file [tmp-dir params name]
            (let [cur-param (get params name)
                  out-file (fs/file tmp-dir (:filename cur-param))]
              [(keyword name) (if (> (:size cur-param) 0)
                                (do
                                  (copy (:tempfile cur-param) out-file)
                                  (str out-file)))]))]
    (let [work-info (prep-tmp-dir)
          in-files (into {} (map (partial download-file (:dir work-info) params)
                                 [:variant-file :region-file]))]
      (session/put! :work-info work-info)
      (session/put! :in-files in-files)
      (scoring-html))))

;; ## File retrieval from processing

(defn get-variant-file
  "Retrieve processed output file for web display."
  [name]
  (letfn [(sample-file [ext]
            (let [base-name "contestant-reference"]
              (format "%s-%s-%s" (get-in @web-config [:ref :sample]) base-name ext)))]
    (let [file-map {"concordant" (sample-file "concordant.vcf")
                    "discordant" (sample-file "discordant.vcf")
                    "discordant-missing" (sample-file "discordant-missing.vcf")
                    "phasing" (sample-file "phasing-error.vcf")}
          base-dir (:dir (session/get :work-info))
          work-dir (if-not (nil? base-dir) (fs/file base-dir "grading"))
          name (get file-map name)
          fname (if-not (or (nil? work-dir)
                            (nil? name)) (str (fs/file work-dir name)))]
      (response/content-type "text/plain"
                             (if (and (not (nil? fname)) (fs/exists? fname))
                               (slurp fname)
                               "Variant file not found")))))
