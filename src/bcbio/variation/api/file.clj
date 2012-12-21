(ns bcbio.variation.api.file
  "Provide top level API for retrieving available files for a user.
  Encapsulates distributed storage in GenomeSpace as well as locally
  produced files."
  (:use [clojure.java.io]
        [bcbio.variation.api.shared :only [web-config remote-file-cache]])
  (:require [clojure.java.shell :as shell]
            [clojure.string :as string]
            [clj-stacktrace.repl :as stacktrace]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.variation.annotate.effects :as effects]
            [bcbio.variation.index.metrics :as metrics]
            [bcbio.variation.index.gemini :as gemini]
            [bcbio.variation.normalize :as normalize]
            [bcbio.variation.remote.core :as remote]))

;; ## File retrieval, with caching

(defn- get-default-public
  "Retrieve default + public directories for specific remote instances."
  [rclient]
  (case (:type rclient)
    :gs (map (fn [x] {:id x})
             (cons "." (remove #(.contains % (:username rclient))
                               (get-in @web-config [:remote :public]))))
    :galaxy (remote/list-dirs rclient nil)
    []))

(defn- update-user-files
  "Update file cache for our current user and filetype"
  [rclient ftype]
   (let [dirnames (get-default-public rclient)
         file-info (mapcat #(remote/list-files rclient % ftype) dirnames)]
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

(def ^{:doc "Set of files currently under preparation, preventing double work."
       :private true}
  prep-queue (atom #{}))

(defn- prep-file*
  "Do actual work of preparing input file: resort and annotate for effects."
  [in-file ref-info out-dir cache-dir]
  (let [current-ref (normalize/pick-best-ref in-file (cons (:genome ref-info) (:genome-alts ref-info)))]
    (-> in-file
        (normalize/prep-vcf (:genome ref-info) nil :out-dir out-dir :orig-ref-file current-ref)
        (effects/snpeff-annotate (:effects ref-info) cache-dir :out-dir out-dir))))

(def ^{:doc "Provide knob to avoid doing labor intensive prep for testing environments."
       :dynamic true}
  *skip-prep* false)

(defn- prep-file
  "Setup preparation for preparing downloading input files."
  [in-file ref-info is-new?]
  (let [out-file (if *skip-prep* in-file (itx/add-file-part in-file "prep"))
        cache-dir (get-in @web-config [:dir :cache])]
    (when (or (itx/needs-run? out-file) is-new?)
      (when-not (contains? @prep-queue in-file)
        (try
          (swap! prep-queue conj in-file)
          (itx/with-temp-dir [out-dir (fs/parent in-file)]
            (fs/rename (prep-file* in-file ref-info out-dir cache-dir)
                       out-file))
          (catch Exception ex
            (stacktrace/pst ex)
            (throw ex))
          (finally
           (swap! prep-queue disj in-file)))))
    out-file))

(def ^{:doc "Provide list of files currently indexing."}
  index-queue (atom #{}))

(defn- is-file-new?
  "Check if a file is newly created, or out of date with server."
  [local-file finfo]
  (if (or (nil? (:created-on finfo)) (nil? (:size finfo)))
    false
    (or (not= (fs/size local-file) (:size finfo))
        (<= (fs/mod-time local-file) (.getTime (:created-on finfo))))))

(defn get-prep-and-index
  "Retrieve a file from remote server and prepare for analysis:
   - Checks for remote updates
   - Downloads the file if necessary
   - Preps the file for analysis, including resorting to reference genome
     and annotating.
   - Indexes for metric retrieval."
  [finfo rclient]
  (let [ref-info (first (:ref @web-config))
        test-local-file (remote/get-file (:id finfo) rclient)
        is-new? (is-file-new? test-local-file finfo)
        local-file (if is-new?
                     (do
                       (fs/delete test-local-file)
                       (remote/get-file (:id finfo) rclient))
                     test-local-file)
        ready-file (prep-file local-file ref-info is-new?)]
    (when ready-file
      (when-not (contains? @index-queue ready-file)
        (try
          (swap! index-queue conj ready-file)
          (metrics/index-variant-file ready-file (:genome ref-info) :re-index? is-new?
                                      :subsample-params (:params @web-config))
          (gemini/index-variant-file ready-file (:genome ref-info) :re-index? is-new?
                                     :subsample-params (:params @web-config))
          (catch Exception ex
            (stacktrace/pst ex))
          (finally
           (swap! index-queue disj ready-file)))))
    ready-file))

(defn pre-fetch-remotes
  "Retrieve and pre-index files for analysis from the remote client."
  [rclient]
  (doall (map (partial update-user-files rclient) [:vcf]))
  (when-let [cache-dir (get-in @web-config [:dir :cache])]
    (when-let [biodata-dir (get-in @web-config [:remote :biodata])]
      (prep-biodata-dir rclient biodata-dir))
    (doseq [x (list-files-w-cache rclient :vcf)]
      (get-prep-and-index x rclient))))

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

