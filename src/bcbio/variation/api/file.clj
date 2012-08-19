(ns bcbio.variation.api.file
  "Provide top level API for retrieving available files for a user.
  Encapsulates distributed storage in GenomeSpace as well as locally
  produced files."
  (:use [clojure.java.io]
        [bcbio.variation.api.shared :only [web-config remote-file-cache]])
  (:require [clojure.string :as string]
            [clj-genomespace.core :as gs]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]))

(declare pre-retrieve-gs)
(defn get-gs-client
  "Top level retrieval of a client from username/password to pre-connected client.
   As a side effect, pre-retrieve and caches files and associated information."
  [creds & {:keys [pre-fetch?]
            :or {pre-fetch? true}}]
  (let [{:keys [username password client]} creds]
    (cond
     (and client (gs/logged-in? client)) (do
                                           (when pre-fetch?
                                             (future (pre-retrieve-gs client)))
                                           client)
     (and username password) (gs/get-client username :password password)
     :else nil)))

(defn- get-gs-dirname-files
  "Retrieve files of the specified type from GenomeSpace."
  [ftype gs-client dirname]
  (concat
   (map (fn [finfo]
          {:id (str "gs:" (:dirname finfo) "/" (:name finfo))
           :tags (remove nil?
                         [(first (drop 3 (string/split (:dirname finfo) #"/" 4)))])
           :folder (:dirname finfo) :filename (:name finfo)
           :created-on (:date finfo)})
        (gs/list-files gs-client dirname (name ftype)))
   (mapcat (partial get-gs-dirname-files ftype gs-client) (gs/list-dirs gs-client dirname))))

(defn get-files
  "Retrieve files information for files on the provided type.
   - creds is a dictionary of login credentials for GenomeSpace and can be:
     :username/:password for login or :client for a connected client.
   - dirnames specify specific directories to list, defaults to all directories."
  [ftype creds & {:keys [dirnames use-cache?]
                  :or {use-cache? true}}]
  (let [dirnames (if (or (nil? dirnames) (empty? dirnames))
                   (cons "." (get-in @web-config [:remote :public]))
                   dirnames)]
    (when-let [gs-client (get-gs-client creds :pre-fetch? use-cache?)]
      (if-let [cache-info (and use-cache? (get @remote-file-cache
                                             [(gs/get-username gs-client) ftype]))]
        cache-info
        (apply concat
               (map (partial get-gs-dirname-files ftype gs-client) dirnames))))))

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
      (let [gs-client (get-gs-client creds)
            out-file (str (file local-dir (fs/base-name remote-name)))]
        (when-not (contains? download-queue out-file)
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
  [fnames base-file subdir creds]
  (let [gs-dir (str (fs/file (fs/parent (last (string/split base-file #":" 2))) subdir))]
    (let [gs-client (get-gs-client creds)]
      (doseq [fname fnames]
        (gs/upload gs-client gs-dir fname)))
    gs-dir))

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
          (doall (map #(retrieve-file (:id %) creds)
                      (get-files :notyet-fa creds :dirnames [biodata-dir]))))))))