(ns bcbio.variation.web.db
  "Provide basic persistence of user files and processes in local DB."
  (:import [com.mchange.v2.c3p0 ComboPooledDataSource])
  (:require [clojure.string :as string]
            [clojure.java.jdbc :as sql]
            [fs.core :as fs]))

(defn get-sqlite-db [fname & {:as opts}]
  "Retrieve SQLite database connection"
  (merge
   {:classname "org.sqlite.JDBC"
    :subprotocol "sqlite"
    :subname fname}
   opts))

(defn get-sqlite-db-pool [fname]
  (let [spec (get-sqlite-db fname)]
    {:datasource (doto (ComboPooledDataSource.)
                   (.setDriverClass (:classname spec))
                   (.setJdbcUrl (str "jdbc:" (:subprotocol spec) ":"
                                     (:subname spec))))}))

(defn- create-user-tables []
  (sql/create-table :analysis
                    [:analysis_id :text "PRIMARY KEY"]
                    [:username :text]
                    [:type :text]
                    [:description :text]
                    [:location :text]
                    [:created :timestamp "NOT NULL" "DEFAULT CURRENT_TIMESTAMP"])
  (sql/create-table :files
                    [:analysis_id :text]
                    [:name :text]
                    [:location :text]))

(defn prepare-web-db
  "Prepare input database for storing user and file information in SQLite."
  [db-file]
  (when-not (fs/exists? (fs/parent db-file))
    (fs/mkdirs (fs/parent db-file)))
  (when-not (fs/exists? db-file)
    (sql/with-connection (get-sqlite-db db-file :create true)
      (sql/transaction
       (create-user-tables))))
  db-file)

(defn get-analyses
  "Retrieve list of analyses run for a specific user and analysis type"
  [username atype db-file]
  (sql/with-connection (get-sqlite-db db-file)
    (sql/with-query-results rows
      ["SELECT * FROM analysis WHERE username = ? AND type = ? ORDER BY created DESC"
       username atype]
      (vec (map #(assoc % :created (java.sql.Timestamp. (:created %))) rows)))))

(defmulti add-analysis
  "Add an analysis and associated files to the database."
  (fn [info _] (:type info)))

(defmethod add-analysis :scoring
  [info db-file]
  (letfn [(get-analysis-files [info]
            (map (fn [[k f]]
                   {:analysis_id (:analysis_id info)
                    :name (name k)
                    :location (string/replace f (str (:location info) "/") "")})
                 (:files info)))]
    (sql/with-connection (get-sqlite-db db-file)
      (sql/transaction
       (sql/insert-record :analysis (-> info
                                        (dissoc :files)
                                        (assoc :created (java.sql.Timestamp. (.getTime (java.util.Date.))))))
       (doseq [x (get-analysis-files info)]
         (sql/insert-record :files x))))))