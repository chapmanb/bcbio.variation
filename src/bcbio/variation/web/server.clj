(ns bcbio.variation.web.server
  "Server providing routes for serving up static pages."
  (:use [clojure.java.io]
        [compojure.core]
        [cemerick.shoreleave.rpc :only [defremote]]
        [ring.adapter.jetty :only [run-jetty]]
        [ring.middleware file file-info keyword-params
         multipart-params nested-params params session]
        [ring.util.response]
        [bcbio.variation.config :only [get-log-status]]
        [bcbio.variation.web.db :only [prepare-web-db get-analyses]]
        [bcbio.variation.web.rpc :only [wrap-rpc-session]]
        [bcbio.variation.api.shared :only [web-config set-config-from-file!]]
        [ring.middleware file anti-forgery file-info])
  (:require [clojure.string :as string]
            [compojure.route :as route]
            [clj-yaml.core :as yaml]
            [fs.core :as fs]
            [bcbio.variation.remote.core :as remote]
            [bcbio.variation.web.process :as web-process]))

(defremote login [creds session]
  (let [rclient (remote/get-client (assoc creds :type :gs))]
    (when (:conn rclient)
      {:result (:username rclient)
       :session (-> session
                    (assoc :rclient rclient))})))

(defremote logout [session]
  {:result nil
   :session (dissoc session :rclient)})

(defn- get-username* [session]
  (get-in session [:rclient :username]))

(defremote get-username [session]
  (get-username* session))

(defremote get-genomes [_]
  (map (fn [x] {:value (:sample x)
                :text (format "%s (%s)" (:sample x) (:description x))})
       (:ref @web-config)))

(defn- prep-display-path [x]
  {:full (:id x)
   :name (last (string/split (:name x) #"/"))})

(defremote list-external-dirs [session]
  (if-let [rclient (:rclient session)]
    (map prep-display-path (remote/list-dirs rclient "."))
    []))

(defremote list-external-files [dir ftype session]
  (if-let [rclient (:rclient session)]
    (map (fn [x] {:full (:id x)
                  :name (:filename x)})
         (remote/list-files rclient {:id dir} ftype))
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
          (future (web-process/run-scoring work-info params (:rclient session)))
          (-> (response out-html)
              (assoc :session
                (assoc session :work-info new-work-info)))))
  (GET "/scorefile/:runid/:name" [runid name :as {session :session}]
       (-> (response (web-process/get-variant-file runid name (get-username* session)
                                                   (get (:work-info session) runid)))
           (content-type "text/plain")))
  (GET "/" req (response (web-process/html-submit-page)))
  (route/files "/" {:root "public" :allow-symlinks? true})
  (route/not-found "Not found"))

(def main-handler
  (-> main-routes
      wrap-file-info
      wrap-anti-forgery
      wrap-rpc-session
      wrap-session
      wrap-keyword-params
      wrap-params
      wrap-multipart-params))

(defn default-config []
  (set-config-from-file! "config/web-processing.yaml"))

(defn -main
  ([config-file]
     (-main config-file "8080"))
  ([config-file port]
     (set-config-from-file! config-file)
     (run-jetty #'main-handler {:join? false :port (Integer/parseInt port)})))
