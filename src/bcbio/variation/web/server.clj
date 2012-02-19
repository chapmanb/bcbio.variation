(ns bcbio.variation.web.server
  "Server providing routes for serving up static pages."
  (:use [clojure.java.io]
        [ring.adapter.jetty :only (run-jetty)]
        [ring.middleware.file :only (wrap-file)]
        [ring.middleware.file-info :only (wrap-file-info)]
        [ring.middleware.params :only (wrap-params)]
        [ring.middleware.multipart-params :only (wrap-multipart-params)]
        [ring.middleware.reload :only (wrap-reload)]
        [ring.middleware.session :only (wrap-session)]
        [ring.util.response :only (file-response redirect)]
        [compojure.core :only (defroutes ANY POST GET)])
  (:require [clj-yaml.core :as yaml]
            [bcbio.variation.web.process :as web-process]))

(def ^:private config (atom nil))

(defroutes app-routes
  (POST "/score" request web-process/prep-scoring)
  (GET "/summary" request web-process/run-scoring)
  (ANY "*" request (file-response "404.html" {:root (-> @config :dir :html-root)})))

(defn wrap-add-config
  "Add configuration information to the current request, loaded from input YAML file."
  [handler]
  (fn [request]
    (handler (assoc request :config @config))))

(defn app []
  (-> app-routes
      (wrap-reload '(web-process))
      (wrap-file (-> @config :dir :html-root))
      wrap-file-info
      wrap-session
      wrap-params
      wrap-multipart-params
      wrap-add-config))

(defn -main
  ([config-file]
     (-main config-file "8080"))
  ([config-file port]
     (reset! config (-> config-file slurp yaml/parse-string))
     (println (str "Running server on http://localhost:" port))
     (run-jetty (app) {:join? false :port (Integer/parseInt port)})))
