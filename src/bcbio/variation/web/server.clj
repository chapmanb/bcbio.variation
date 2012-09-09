(ns bcbio.variation.web.server
  "Server providing routes for serving up static pages."
  (:use [clojure.java.io]
        [compojure.core]
        [cemerick.shoreleave.rpc :only [defremote]]
        [ring.adapter.jetty :only [run-jetty]]
        [ring.middleware file file-info keyword-params
         multipart-params nested-params params session]
        [ring.util.response]
        [bcbio.variation.api.file :only [get-files]]
        [bcbio.variation.config :only [get-log-status]]
        [bcbio.variation.web.db :only [prepare-web-db get-analyses]]
        [bcbio.variation.web.rpc :only [wrap-rpc-session]]
        [bcbio.variation.api.shared :only [web-config set-config-from-file!]]
        [ring.middleware file anti-forgery file-info])
  (:require [clojure.string :as string]
            [compojure.route :as route]
            [clj-yaml.core :as yaml]
            [fs.core :as fs]
            [clj-genomespace.core :as gs]
            [bcbio.variation.web.process :as web-process]))

(defn- cur-gs-client [session]
  (when-let [gs-client (:gs-client session)]
    (when (gs/logged-in? gs-client)
      gs-client)))

(defremote login [{:keys [username password]} session]
  (when-let [gs-client (gs/get-client username :password password)]
    {:result username
     :session (-> session
                  (assoc :username username)
                  (assoc :gs-client gs-client))}))

(defremote logout [session]
  {:result nil
   :session (-> session (dissoc :username) (dissoc :gs-client))})

(defn- get-username* [session]
  (when (cur-gs-client session)
    (:username session)))

(defremote get-username [session]
  (get-username* session))

(defremote get-genomes [_]
  (map (fn [x] {:value (:sample x)
                :text (format "%s (%s)" (:sample x) (:description x))})
       (:ref @web-config)))

(defn- prep-gs-path [x]
  {:full x
   :name (last (string/split x #"/"))})

(defremote list-external-dirs [session]
  (if-let [gs-client (cur-gs-client session)]
    (map prep-gs-path (gs/list-dirs gs-client "."))
    []))

(defremote list-external-files [dir ftype session]
  (if-let [gs-client (cur-gs-client session)]
    (map (fn [x] {:full (str (:folder x) "/" (:filename x))
                  :name (:filename x)})
         (get-files ftype {:client gs-client} :dirnames [dir]
                    :use-cache? false))
    []))

(defremote get-status [run-id session]
  (get-log-status {:dir {:out (-> (:work-info session)
                                  (get run-id)
                                  :dir
                                  (file "grading"))}}))

(defremote get-summary [run-id session]
  (when-let [out-dir (if-let [username (get-username* session)]
                       (->> (get-analyses username :scoring (:db @web-config))
                            (filter #(= run-id (:analysis_id %)))
                            first
                            :location)
                       (-> (:work-info session)
                           (get run-id)
                           :dir))]
    (let [summary-file (file out-dir "scoring-summary.html")]
      (when (fs/exists? summary-file)
        (slurp summary-file)))))

(defroutes main-routes
  (GET "/analyses" [:as {session :session}]
       (response (web-process/analyses-html (get-username* session))))
  (POST "/score" [:as {params :params session :session}]
        (let [{:keys [work-info out-html]} (web-process/prep-scoring params)
              new-work-info (assoc (:work-info session)
                              (:id work-info) work-info)]
          (future (web-process/run-scoring work-info params
                                           (get-username* session)
                                           (cur-gs-client session)))
          (-> (response out-html)
              (assoc :session
                (assoc session :work-info new-work-info)))))
  (GET "/scorefile/:runid/:name" [runid name :as {session :session}]
       (-> (response (web-process/get-variant-file runid name (get-username* session)
                                                   (get (:work-info session) runid)))
           (content-type "text/plain")))
  (route/files "/" {:root "public" :allow-symlinks? true})
  (route/not-found "Not found"))

(defn app []
  (-> main-routes
      (wrap-file (get-in @web-config [:dir :html-root]))
      wrap-rpc-session
      wrap-session
      wrap-keyword-params
      wrap-params
      wrap-multipart-params
      ;wrap-file-info
      ;wrap-anti-forgery
      ))

(defn -main
  ([config-file]
     (-main config-file "8080"))
  ([config-file port]
     (set-config-from-file! config-file)
     (run-jetty (app) {:join? false :port (Integer/parseInt port)})))
