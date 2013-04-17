(ns bcbio.variation.api.shared
  "Shared functionality useful across multiple API calls."
  (:use [clojure.java.io]
        [bcbio.variation.remote.client :only [gs-default-server]]
        [bcbio.variation.web.db :only [prepare-web-db]])
  (:require [clojure.string :as string]
            [clj-yaml.core :as yaml]))

(def ^{:doc "Web configuration, loaded from input YAML file"}
  web-config (atom nil))

(def ^{:doc "Hold directory of remote files by user and filetype."}
  remote-file-cache (atom {}))

(defn url->dir
  [url]
  (string/replace (.getHost (as-url url)) "." "_"))

(defn load-web-config
  [config-file]
  (let [config (-> config-file slurp yaml/parse-string)]
    (letfn [(maybe-fix-biodata [x]
              (if (.startsWith x "biodata:")
                (str (get-in config [:dir :cache])
                     "/" (url->dir gs-default-server)
                     (get-in config [:remote :biodata])
                     (string/replace-first x "biodata:" ""))
                x))
            (fix-gs-ref [ref]
              (reduce (fn [coll k]
                        (let [val (get coll k)]
                          (assoc coll k
                                 (cond
                                  (string? val) (maybe-fix-biodata val)
                                  (seq? val) (map maybe-fix-biodata val)
                                  :else val))))
                      ref (keys ref)))]
      (assoc config :ref (map fix-gs-ref (:ref config))))))

(defn set-config-from-file!
  "Set configuration and database information from input YAML file."
  [config-file]
  (let [config (load-web-config config-file)]
    (reset! web-config (assoc config :db
                              (prepare-web-db (str (file (get-in config [:dir :work])
                                                         "analyses.db")))))))
