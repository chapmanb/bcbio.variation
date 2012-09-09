(ns bcbio.variation.web.rpc
  "Provide shoreleave compatible RPC for ring, handling sessions.
  Modified from Chas Emerick: https://github.com/cemerick/shoreleave-remote-ring"
  (require [cemerick.shoreleave.rpc :as baserpc]))

(defn- map-w-session? [x]
  (and (instance? clojure.lang.PersistentArrayMap x)
       (contains? x :session)))

(defn call-remote-session
  [remote-key params session]
  (if-let [func (@baserpc/remotes remote-key)]
    (let [out (apply func (conj (vec params) session))
          result (if (map-w-session? out) (:result out) out)
          session (if (map-w-session? out) (:session out) session)]
      {:status 202
       :headers {"Content-Type" "application/clojure; charset=utf-8"}
       :session session
       :body (pr-str result)})
    {:status 404}))

(defn handle-rpc
  [{{:keys [params remote]} :params :as request}]
  (call-remote-session (keyword remote) (baserpc/safe-read params)
                       (:session request)))

(defn wrap-rpc-session
  ([app] (wrap-rpc-session app baserpc/default-remote-uri))
  ([app remote-uri]
     (fn [{:keys [request-method uri] :as request}]
       (if (and (= :post request-method) (= remote-uri uri))
         (handle-rpc request)
         (app request)))))
