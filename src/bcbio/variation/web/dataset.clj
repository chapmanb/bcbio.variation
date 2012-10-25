(ns bcbio.variation.web.dataset
  "Provide retrieval of local file datasets, for display and upload to other services."
  (:import [java.net InetAddress URL]
           [java.util UUID])
  (:use [clojure.java.io]
        [bcbio.variation.api.shared :only [web-config]])
  (:require [clojure.string :as string]
            [fs.core :as fs]
            [bcbio.variation.web.db :as db]))

;; ## Remote file access

(def ^{:private true
       :doc "List of available datasets for external retrieval."}
  exposed-datasets (atom {}))

(defn expose
  "Provide a remote file for a single download via a remote server.
   Returns dataset identifier to use for retrieval."
  [fname remote-url]
  {:pre [(fs/exists? fname)]}
  (let [dsid (str (UUID/randomUUID))
        remote-host (.getHost (URL. remote-url))
        expected-remote (set (map #(.getHostAddress %) (InetAddress/getAllByName remote-host)))]
    (swap! exposed-datasets assoc dsid {:fname fname :expected-remote expected-remote})
    dsid))

(defn expose-w-url
  "Expose a dataset with a provided callback URL."
  [fname remote-url cb-host cb-port cb-path]
  (let [dsid (expose fname remote-url)]
    (str "http://" cb-host ":" cb-port "/" cb-path "/" dsid)))

(defn retrieve
  "Retrieve a dataset via identifier, checking remote host for permissions match."
  [dsid remote-addr]
  (when-let [{:keys [fname expected-remote]} (get exposed-datasets dsid)]
    (when (contains? expected-remote remote-addr)
      (swap! exposed-datasets dissoc dsid)
      {:status 200
       :header {}
       :body (input-stream fname)})))

;; ## File retrieval from processing

(defn- get-run-info
  "Retrieve run information from stored database or current work-info."
  [run-id username session-work-info]
  (if (nil? username)
    [(:comparison-genome session-work-info) (:dir session-work-info)]
    (let [work-info (->> (db/get-analyses username :scoring (:db @web-config))
                         (filter #(= run-id (:analysis_id %)))
                         first)
          sample-name (when-not (nil? work-info)
                        (first (string/split (:description work-info) #":")))]
      [sample-name (:location work-info)])))

(defn get-variant-file
  "Retrieve processed output file for web display."
  [run-id name username work-info]
  (letfn [(sample-file [sample-name ext]
            (let [base-name "contestant-reference"]
              (format "%s-%s-%s" sample-name base-name ext)))]
    (let [[sample-name base-dir] (get-run-info run-id username work-info)
          file-map {"concordant" (sample-file sample-name "concordant.vcf")
                    "discordant" (sample-file sample-name "discordant.vcf")
                    "discordant-missing" (sample-file sample-name "discordant-missing.vcf")
                    "phasing" (sample-file sample-name "phasing-error.vcf")}
          work-dir (when-not (nil? base-dir) (fs/file base-dir "grading"))
          name (get file-map name)
          fname (if-not (or (nil? work-dir)
                            (nil? name)) (str (fs/file work-dir name)))]

      (if (and (not (nil? fname)) (fs/exists? fname))
        (slurp fname)
        "Variant file not found"))))