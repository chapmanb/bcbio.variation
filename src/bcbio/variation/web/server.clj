(ns bcbio.variation.web.server
  "Server providing routes for serving up static pages."
  (:use [clojure.java.io]
        [noir.core :only [defpage]]
        [noir.fetch.remotes :only [defremote]]
        [bcbio.variation.api.file :only [get-files]]
        [bcbio.variation.config :only [get-log-status]]
        [bcbio.variation.web.db :only [prepare-web-db get-analyses]]
        [bcbio.variation.api.shared :only [web-config set-config-from-file!]]
        [ring.middleware file anti-forgery file-info])
  (:require [clojure.string :as string]
            [clj-yaml.core :as yaml]
            [fs.core :as fs]
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
    (map (fn [x] {:full (str (:folder x) "/" (:filename x))
                  :name (:filename x)})
         (get-files ftype {:client gs-client} [dir]))
    []))

(defremote get-status [run-id]
  (get-log-status {:dir {:out (-> (session/get :work-info)
                                  (get run-id)
                                  :dir
                                  (file "grading"))}}))

(defremote get-summary [run-id]
  (when-let [out-dir (if-let [username (get-username)]
                       (->> (get-analyses username :scoring (:db @web-config))
                            (filter #(= run-id (:analysis_id %)))
                            first
                            :location)
                       (-> (session/get :work-info)
                           (get run-id)
                           :dir))]
    (let [summary-file (file out-dir "scoring-summary.html")]
      (when (fs/exists? summary-file)
        (slurp summary-file)))))

(defpage "/analyses" {}
  (web-process/analyses-html (get-username)))

(defpage [:post "/score"] {:as params}
  (let [{:keys [work-info out-html]} (web-process/prep-scoring params)]
    (future (web-process/run-scoring work-info params))
    out-html))

(defpage "/scorefile/:runid/:name" {:keys [runid name]}
  (web-process/get-variant-file runid name (get-username)))

(defn -main
  ([config-file]
     (-main config-file "8080"))
  ([config-file port]
     (set-config-from-file! config-file)
     (server/add-middleware wrap-file (get-in @web-config [:dir :html-root]))
     ;;(server/add-middleware wrap-file-info)
     ;;(server/add-middleware wrap-anti-forgery)
     (server/start (Integer/parseInt port))))
