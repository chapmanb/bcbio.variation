(ns bcbio.variation.api.file
  "Provide top level API for retrieving available files for a user.
  Encapsulates distributed storage in GenomeSpace as well as locally
  produced files."
  (:use [clojure.java.io])
  (:require [clojure.string :as string]
            [clj-genomespace.core :as gs]
            [fs.core :as fs]))

(defn get-gs-client
  "Top level retrieval of a client from username/password to pre-connected client."
  [creds]
  (let [{:keys [username password client]} creds]
    (cond
     (and client (gs/logged-in? client)) client
     (and username password) (gs/get-client username :password password)
     :else nil)))

(defn- get-gs-dirname-files
  "Retrieve files of the specified type from GenomeSpace."
  [ftype gs-client dirname]
  (map (fn [finfo]
         {:id (str "gs:" (:dirname finfo) "/" (:name finfo))
          :tags (remove nil?
                        [(first (drop 3 (string/split (:dirname finfo) #"/" 4)))])
          :folder (:dirname finfo) :filename (:name finfo)
          :created-on (:date finfo)})
       (gs/list-files gs-client dirname (name ftype))))

(defn get-files
  "Retrieve files information for files on the provided type.
   - creds is a dictionary of login credentials for GenomeSpace and can be:
     :username/:password for login or :client for a connected client.
   - dirnames specify specific directories to list, defaults to all directories."
  ([ftype creds]
     (get-files ftype creds nil))
  ([ftype creds dirnames]
     (when-let [gs-client (get-gs-client creds)]
       (let [dirnames (if (or (nil? dirnames) (empty? dirnames))
                        (cons "." (gs/list-dirs gs-client "."))
                        dirnames)]
         (apply concat
                (map (partial get-gs-dirname-files ftype gs-client) dirnames))))))

(defmulti retrieve-file
  "Retrieve files by name, transparently handling remote files."
  (fn [fname _ _]
    (let [parts (string/split fname #":" 2)]
      (when (= 2 (count parts))
        (keyword (first parts))))))

(defmethod retrieve-file :gs
  [fname creds cache-dir]
  (let [remote-name (second (string/split fname #":" 2))
        local-file (str (file cache-dir (if (.startsWith remote-name "/")
                                          (subs remote-name 1)
                                          remote-name)))
        local-dir (str (fs/parent local-file))]
    (when-not (fs/exists? local-file)
      (when-not (fs/exists? local-dir)
        (fs/mkdirs local-dir))
      (let [gs-client (get-gs-client creds)]
        (gs/download gs-client (str (fs/parent remote-name))
                     (fs/base-name remote-name) local-dir)))
    local-file))

(defmethod retrieve-file :default [fname _ _]
  fname)

(defn put-files
  "Put files back on remote GenomeSpace server relative to the base file."
  [fnames base-file subdir creds]
  (let [gs-dir (str (fs/file (fs/parent (last (string/split base-file #":" 2))) subdir))]
    (let [gs-client (get-gs-client creds)]
      (doseq [fname fnames]
        (gs/upload gs-client gs-dir fname)))
    gs-dir))