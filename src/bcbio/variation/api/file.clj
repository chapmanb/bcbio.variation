(ns bcbio.variation.api.file
  "Provide top level API for retrieving available files for a user.
  Encapsulates distributed storage in GenomeSpace as well as locally
  produced files."
  (:use [clojure.java.io]
        [bcbio.variation.api.shared :only [web-config remote-file-cache]])
  (:require [clojure.java.shell :as shell]
            [clojure.string :as string]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.variation.index.metrics :as metrics]
            [bcbio.variation.index.gemini :as gemini]
            [bcbio.variation.remote.core :as remote]))

;; ## File retrieval, with caching

(defn- update-user-files
  "Update file cache for our current user and filetype"
  [rclient ftype]
   (let [dirnames (cons "." (remove #(.contains % (:username rclient))
                                    (get-in @web-config [:remote :public])))
         file-info (mapcat #(remote/list-files rclient {:id %} ftype) dirnames)]
     (swap! remote-file-cache assoc [(:username rclient) ftype] file-info)
     file-info))

(defn list-files-w-cache
  "Retrieve file information for files of the specified type, with caching."
  [rclient ftype]
  (let [cache-info (get @remote-file-cache [(:username rclient) ftype])]
    (if (seq cache-info)
      (do
        (future (update-user-files rclient ftype))
        cache-info)
      (update-user-files rclient ftype))))

;; ## Pre-fetching of data for front end interactivity

(defn- prep-biodata-dir
  "Prepare biodata directory, synchronizing with GenomeSpace client."
  [rclient biodata-dir]
  (letfn [(download-and-unpack [gs-id]
            (let [zip-file (remote/get-file gs-id rclient)]
              (when (and (fs/exists? zip-file)
                         (itx/needs-run? (itx/remove-zip-ext zip-file)))
                (shell/sh "gunzip" zip-file)
                (spit zip-file (str "unzipped to " (itx/remove-zip-ext zip-file))))))]
    (doall (map #(download-and-unpack (:id %))
                (remote/list-files rclient biodata-dir :gz)))))

(defn- retrieve-file-check-update
  "Retrieve a file from GenomeSpace checking to see if the file changed remotely."
  [rclient finfo]
  (let [local-file (remote/get-file (:id finfo) rclient)
        is-current? (and (= (fs/size local-file) (:size finfo))
                         (or (nil? (:created-on finfo))
                             (>= (fs/mod-time local-file) (.getTime (:created-on finfo)))))]
    (when-not is-current?
      (fs/delete local-file)
      (remote/get-file (:id finfo) rclient))
    [local-file (not is-current?)]))

(def ^{:doc "Provide list of files currently indexing."}
  index-queue (atom #{}))

(defn pre-index-variants
  "Provide download and pre-indexing of GenomeSpace variant files."
  [rclient index-type index-fn]
  (letfn [(do-pre-index [finfo ref-file]
            (let [[local-file is-new?] (retrieve-file-check-update rclient finfo)
                  k {:type index-type :file local-file}]
              (when-not (contains? @index-queue k)
                (try
                  (swap! index-queue conj k)
                  (index-fn local-file ref-file :re-index? is-new?
                            :subsample-params (:params @web-config))
                  (finally
                   (swap! index-queue disj k))))))]
    (let [ref-file (:genome (first (:ref @web-config)))]
      (doall (map #(do-pre-index % ref-file)
                  (list-files-w-cache rclient :vcf))))))

(defn pre-fetch-remotes
  "Retrieve and pre-index files for analysis from the remote client."
  [rclient]
  (doall (map (partial update-user-files rclient) [:vcf]))
  (when-let [cache-dir (get-in @web-config [:dir :cache])]
    (when-let [biodata-dir (get-in @web-config [:remote :biodata])]
      (prep-biodata-dir rclient biodata-dir))
    (pre-index-variants rclient "metrics" metrics/index-variant-file)
    (pre-index-variants rclient "gemini" gemini/index-variant-file)))

;; ## Client API, with pre-fetching

(defn get-client
  "Top level retrieval of a client from username/password to pre-connected client.
   As a side effect, pre-retrieve and caches files and associated information."
  [creds & {:keys [pre-fetch? allow-offline?]
            :or {pre-fetch? true}}]
  (let [rclient (remote/get-client (-> creds
                                       (assoc :allow-offline? allow-offline?)
                                       (assoc :type (get creds :type :gs))))]
    (when (and pre-fetch? (:conn rclient))
      (future (pre-fetch-remotes rclient)))
    rclient))

(defn- get-local-dl-files
  "Retrieve local cached files supporting offline processing
   XXX Needs updating to hook back in for full offline analysis reboot."
  [ftype & {:keys [dirnames]}]
  {:pre [(nil? dirnames)]}
  (letfn [(check-dir-for-type [root _ files]
            (->> files
                 (filter #(.endsWith % (str "." (name ftype))))
                 (map #(str (file root %)))))
          (convert-to-api [cache-dir fname]
            {:id fname
             :tags []
             :filename (str (fs/base-name fname))
             :folder (string/replace (str (fs/parent fname)) (str (fs/file cache-dir)) "")
             :size (fs/size fname)
             :created-on (java.util.Date. (fs/mod-time fname))})]
    (let [cache-dir (get-in @web-config [:dir :cache])]
      (->> (fs/walk check-dir-for-type cache-dir)
           flatten
           (map (partial convert-to-api cache-dir))))))

