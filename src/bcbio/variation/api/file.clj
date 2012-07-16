(ns bcbio.variation.api.file
  "Provide top level API for retrieving available files for a user.
  Encapsulates distributed storage in GenomeSpace as well as locally
  produced files."
  (require [clj-genomespace.core :as gs]))

(defn- get-gs-client
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
         {:id (str (:dirname finfo) "/" (:name finfo))
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