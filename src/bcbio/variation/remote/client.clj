(ns bcbio.variation.remote.client
  "Establish a connection with a remote service for file management."
  (:require [blend.galaxy.core :as galaxy]
            [clj-genomespace.core :as gs]))

(defrecord RemoteClient [type conn username])

(defmulti get-client
  "Retrieve a remote client using authentication credentials
   creds is a map containing the :type of connection to establish
   and client specific authentication details."
  (fn [creds]
    (:type creds)))

(defmethod get-client :gs
  ^{:doc "Retrieve a GenomeSpace client connection"}
  [creds]
  (let [{:keys [username password client allow-offline?]} creds
        gs-client (cond
                   (and client (gs/logged-in? client)) client
                   (and username password) (try (gs/get-client username :password password)
                                                (catch Exception e
                                                  (when-not allow-offline?
                                                    (throw e))))
                   :else nil)
        username (when gs-client
                   (gs/get-username gs-client))]
    (RemoteClient. :gs gs-client username)))

(defmethod get-client :galaxy
  ^{:doc "Retrieve a Galaxy client connection."}
  [creds]
  (let [{:keys [url api-key client allow-offline?]} creds
        galaxy-client (cond
                       (not (nil? client)) client
                       (and url api-key) (galaxy/get-client url api-key)
                       :else nil)
        user-info (try (galaxy/get-user-info galaxy-client)
                       (catch Exception e
                         (when-not allow-offline?
                           (throw e))))
        username (when user-info
                   (get user-info :username (:email user-info)))]
    (RemoteClient. :galaxy (when user-info galaxy-client) username)))