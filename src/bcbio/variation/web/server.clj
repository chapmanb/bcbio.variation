(ns bcbio.variation.web.server
  "Server providing routes for serving up static pages."
  (:use [clojure.java.io]
        [noir.core :only [defpage]]
        [noir.fetch.remotes :only [defremote]]
        [bcbio.variation.web.shared :only [web-config]]
        [ring.middleware file])
  (:require [clojure.string :as string]
            [clj-yaml.core :as yaml]
            [noir.server :as server]
            [noir.session :as session]
            [clj-genomespace.core :as gs]
            [bcbio.variation.web.process :as web-process]))

(defn- cur-gs-client []
  (when-let [gs-client (session/get :gs-client)]
    (when (gs/logged-in? gs-client)
      gs-client)))

(defremote login [{:keys [username password]}]
  (when-let [gs-client (gs/get-client username :password password)]
    (session/put! :username username)
    (session/put! :gs-client gs-client)
    username))

(defremote logout []
  (session/remove! :username)
  (session/remove! :gs-client)
  nil)

(defremote get-username []
  (when (cur-gs-client)
    (session/get :username)))

(defremote get-genomes []
  (map (fn [x] {:value (:sample x)
                :text (format "%s (%s)" (:sample x) (:description x))})
       (:ref @web-config)))

(defn- prep-gs-path [x]
  {:full x
   :name (last (string/split x #"/"))})

(defremote list-external-dirs []
  (if-let [gs-client (cur-gs-client)]
    (map prep-gs-path (gs/list-dirs gs-client "."))
    []))

(defremote list-external-files [dir ftype]
  (if-let [gs-client (cur-gs-client)]
    (map prep-gs-path (gs/list-files gs-client dir ftype))
    []))

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
