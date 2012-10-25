(ns bcbio.variation.web.dataset
  "Provide retrieval of local file datasets, for display and upload to other services."
  (:import [java.net InetAddress URL]
           [java.util UUID])
  (:use [clojure.java.io])
  (:require [fs.core :as fs]))

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