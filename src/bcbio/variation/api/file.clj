(ns bcbio.variation.api.file
  "Provide top level API for retrieving available files for a user.
  Encapsulates distributed storage in GenomeSpace as well as locally
  produced files."
  (:use [clojure.java.io]
        [bcbio.variation.api.shared :only [web-config remote-file-cache]])
  (:require [clojure.java.shell :as shell]
            [clojure.string :as string]
            [clj-genomespace.core :as gs]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.variation.index.metrics :as metrics]
            [bcbio.variation.index.gemini :as gemini]))

(declare pre-retrieve-gs)
(defn get-gs-client
  "Top level retrieval of a client from username/password to pre-connected client.
   As a side effect, pre-retrieve and caches files and associated information."
  [creds & {:keys [pre-fetch? allow-offline?]
            :or {pre-fetch? true}}]
  (let [{:keys [username password client]} creds
        gs-client (cond
                   (and client (gs/logged-in? client)) client
                   (and username password) (try (gs/get-client username :password password)
                                                (catch Exception e
                                                  (when-not allow-offline?
                                                    (throw e))))
                   :else nil)]
    (when (and pre-fetch? gs-client)
      (future (pre-retrieve-gs client)))
    gs-client))

(defn- get-gs-dirname-files
  "Retrieve files of the specified type from GenomeSpace."
  [ftype gs-client dirname]
  (concat
   (map (fn [finfo]
          {:id (str "gs:" (:dirname finfo) "/" (:name finfo))
           :tags (remove nil?
                         [(first (drop 3 (string/split (:dirname finfo) #"/" 4)))])
           :folder (:dirname finfo) :filename (:name finfo)
           :size (:size finfo)
           :created-on (:date finfo)})
        (gs/list-files gs-client dirname (name ftype)))
   (mapcat (partial get-gs-dirname-files ftype gs-client) (gs/list-dirs gs-client dirname))))

(defn- get-local-dl-files
  "Retrieve local files downloaded from remote GenomeSpace."
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

(defn get-files
  "Retrieve files information for files on the provided type.
   - creds is a dictionary of login credentials for GenomeSpace and can be:
     :username/:password for login or :client for a connected client.
   - dirnames specify specific directories to list, defaults to all directories."
  [ftype creds & {:keys [dirnames use-cache?]
                  :or {use-cache? true}}]
  (if-let [gs-client (get-gs-client creds :pre-fetch? use-cache?)]
    (if-let [cache-info (and use-cache? (nil? dirnames)
                             (get @remote-file-cache
                                  [(gs/get-username gs-client) ftype]))]
      cache-info
      (let [dirnames (if (or (nil? dirnames) (empty? dirnames))
                       (cons "." (remove #(.contains % (gs/get-username gs-client))
                                         (get-in @web-config [:remote :public])))
                       dirnames)]
        (apply concat
               (map (partial get-gs-dirname-files ftype gs-client) dirnames))))
    (get-local-dl-files ftype :dirnames dirnames)))

(def ^{:doc "Provide list of files currently under download."}
  download-queue (atom #{}))

(defmulti retrieve-file
  "Retrieve files by name, transparently handling remote files."
  (fn [fname _]
    (let [parts (string/split fname #":" 2)]
      (when (= 2 (count parts))
        (keyword (first parts))))))

(defmethod retrieve-file :gs
  [fname creds]
  (let [cache-dir (get-in @web-config [:dir :cache])
        remote-name (second (string/split fname #":" 2))
        local-file (str (file cache-dir (if (.startsWith remote-name "/")
                                          (subs remote-name 1)
                                          remote-name)))
        local-dir (str (fs/parent local-file))]
    (when-not (fs/exists? local-file)
      (when-not (fs/exists? local-dir)
        (fs/mkdirs local-dir))
      (let [gs-client (get-gs-client creds :pre-fetch? false)
            out-file (str (file local-dir (fs/base-name remote-name)))]
        (when-not (contains? @download-queue out-file)
          (itx/with-tx-file [out-file-tx out-file]
            (swap! download-queue conj out-file)
            (try
              (gs/download gs-client (str (fs/parent remote-name))
                           (fs/base-name remote-name) out-file-tx)
              (finally
               (swap! download-queue disj out-file)))))))
    local-file))

(defmethod retrieve-file :default [fname _]
  fname)

(defn put-files
  "Put files back on remote GenomeSpace server relative to the base file."
  [fnames base-file subdir creds & {:keys [pre-fetch?]
                                    :or {pre-fetch? true}}]
  (let [gs-dir (str (fs/file (fs/parent (last (string/split base-file #":" 2))) subdir))]
    (let [gs-client (get-gs-client creds :pre-fetch? pre-fetch?)]
      (doseq [fname fnames]
        (gs/upload gs-client gs-dir fname)))
    gs-dir))

;; ## Pre-fetching of data for front end interactivity

(defn- prep-biodata-dir
  "Prepare biodata directory, synchronizing with GenomeSpace client."
  [creds biodata-dir]
  (letfn [(download-and-unpack [gs-id creds]
            (let [zip-file (retrieve-file gs-id creds)]
              (when (and (fs/exists? zip-file)
                         (itx/needs-run? (itx/remove-zip-ext zip-file)))
                (shell/sh "gunzip" zip-file)
                (spit zip-file (str "unzipped to " (itx/remove-zip-ext zip-file))))))]
  (doall (map #(download-and-unpack (:id %) creds)
              (get-files :gz creds :dirnames [biodata-dir]
                         :use-cache? false)))))

(defn- retrieve-file-check-update
  "Retrieve a file from GenomeSpace checking to see if the file changed remotely."
  [finfo creds]
  (let [local-file (retrieve-file (:id finfo) creds)
        is-current? (and (= (fs/size local-file) (:size finfo))
                         (>= (fs/mod-time local-file) (.getTime (:created-on finfo))))]
    (when-not is-current?
      (fs/delete local-file)
      (retrieve-file (:id finfo) creds))
    [local-file (not is-current?)]))

(def ^{:doc "Provide list of files currently indexing."}
  index-queue (atom #{}))

(defn pre-index-variants
  "Provide download and pre-indexing of GenomeSpace variant files."
  [creds index-type index-fn]
  (letfn [(do-pre-index [finfo ref-file]
            (let [[local-file is-new?] (retrieve-file-check-update finfo creds)
                  k {:type index-type :file local-file}]
              (when-not (contains? @index-queue k)
                (try
                  (swap! index-queue conj k)
                  (index-fn local-file ref-file :re-index? is-new?
                            :subsample-params (get-in @web-config [:params]))
                  (finally
                   (swap! index-queue disj k))))))]
    (let [ref-file (:genome (first (:ref @web-config)))]
      (doall (map #(do-pre-index % ref-file)
                  (get-files :vcf creds :use-cache? false))))))

(defn pre-retrieve-gs
  "Retrieve and pre-index files for analysis from GenomeSpace."
  [client]
  (letfn [(update-user-files [ftype creds]
            (swap! remote-file-cache assoc [(gs/get-username client) ftype]
                   (get-files ftype creds :use-cache? false)))]
    (let [creds {:client client}]
      (doall (map #(update-user-files % creds) [:vcf]))
      (when-let [cache-dir (get-in @web-config [:dir :cache])]
        (when-let [biodata-dir (get-in @web-config [:remote :biodata])]
          (prep-biodata-dir creds biodata-dir))
        (pre-index-variants creds "metrics" metrics/index-variant-file)
        (pre-index-variants creds "gemini" gemini/index-variant-file)))))