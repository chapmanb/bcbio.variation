(ns bcbio.variation.api.shared
  "Shared functionality useful across multiple API calls."
  (:use [clojure.java.io]
        [bcbio.variation.web.db :only [prepare-web-db]])
  (:require [clj-yaml.core :as yaml]))
  
(def ^{:doc "Web configuration, loaded from input YAML file"}
  web-config (atom nil))

(defn set-config-from-file!
  "Set configuration and database information from input YAML file."
  [config-file]
  (let [config (-> config-file slurp yaml/parse-string)]
    (reset! web-config (assoc config :db
                              (prepare-web-db (str (file (get-in config [:dir :work])
                                                         "analyses.db")))))))
