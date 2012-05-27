(ns bcbio.variation.web.server
  "Server providing routes for serving up static pages."
  (:use [clojure.java.io]
        [noir.core :only [defpage]]
        [bcbio.variation.web.shared :only [web-config]]
        [ring.middleware file])
  (:require [clj-yaml.core :as yaml]
            [noir.server :as server]
            [bcbio.variation.web.process :as web-process]))

(def ^:private test-usernames
  {"tester" "tester"})

(defpage [:post "/login"] {:keys [username password]}
  (when (= (get test-usernames username) password)
    username))

(defpage [:post "/score"] {:as params}
  (web-process/prep-scoring params))

(defpage "/summary" []
  (web-process/run-scoring))

(defpage "/scorefile/:name" {:keys [name]}
  (web-process/get-variant-file name))

(defn -main
  ([config-file]
     (-main config-file "8080"))
  ([config-file port]
     (reset! web-config (-> config-file slurp yaml/parse-string))
     (server/add-middleware wrap-file (get-in @web-config [:dir :html-root]))
     (server/start (Integer/parseInt port))))
