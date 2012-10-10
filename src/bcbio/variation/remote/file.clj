(ns bcbio.variation.remote.file
  "List, retrieve and push files from a remote filestore."
  (:use [clojure.java.io]
        [bcbio.variation.api.shared :only [web-config]])
  (:require [clojure.string :as string]
            [fs.core :as fs]
            [blend.galaxy.core :as galaxy]
            [clj-genomespace.core :as gs]
            [bcbio.run.itx :as itx]))

(def ^{:doc "Provide list of files currently under download."
       :private true}
  download-queue (atom #{}))

(defn- download-to-local
  "Generalized remote download to local cache directory.
   download-fn takes the client, remote name and local name, and
   handles the download step for a specific remote service."
  [fname rclient download-fn]
  (let [cache-dir (get-in @web-config [:dir :cache])
        remote-name (second (string/split fname #":" 2))
        local-file (str (file cache-dir (if (.startsWith remote-name "/")
                                          (subs remote-name 1)
                                          remote-name)))
        local-dir (str (fs/parent local-file))]
    (when-not (fs/exists? local-file)
      (when-not (fs/exists? local-dir)
        (fs/mkdirs local-dir))
      (let [out-file (str (file local-dir (fs/base-name remote-name)))]
        (when-not (contains? @download-queue out-file)
          (itx/with-tx-file [out-file-tx out-file]
            (swap! download-queue conj out-file)
            (try
              (download-fn rclient remote-name out-file-tx)
              (finally
               (swap! download-queue disj out-file)))))))
    local-file))

(defmulti get-file
  "Retrieve files by name, transparently handling remote files."
  (fn [fname _]
    (let [parts (string/split fname #":" 2)]
      (when (= 2 (count parts))
        (keyword (first parts))))))

(defmethod get-file :gs
  ^{:doc "Retrieve a file from GenomeSpace to the local cache"}
  [fname rclient]
  (letfn [(download-gs [rclient remote-name out-file]
            (gs/download (:client rclient) (str (fs/parent remote-name))
                         (fs/base-name remote-name) out-file))]
    (download-to-local fname rclient download-gs)))

(defmethod get-file :galaxy
  ^{:doc "Retrieve a file from Galaxy to the local cache"}
  [fname rclient]
  (letfn [(download-galaxy [rclient remote-name out-file]
            (let [[_ history-id fname] (string/split remote-name #"/")]))]
    (download-to-local fname rclient download-galaxy)))

(defmethod get-file :default
  ^{:doc "Get local file: no-op, just return the file."}
  [fname _]
  fname)
