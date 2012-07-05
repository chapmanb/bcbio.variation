(ns bcbio.variation.web.db
  "Provide basic persistence of user files and processes in local DB."
  (:require [clojure.java.jdbc :as sql]
            [fs.core :as fs]))

(defn- get-db [fname & {:as opts}]
  "Retrieve SQLite database connection"
  (merge
   {:classname "org.sqlite.JDBC"
    :subprotocol "sqlite"
    :subname fname}
   opts))

(defn- create-user-tables []
  (sql/create-table :analysis
                    [:analysis_id :text "PRIMARY KEY"]
                    [:username :text]
                    [:type :text]
                    [:description :text]
                    [:created :timestamp "NOT NULL" "DEFAULT CURRENT_TIMESTAMP"])
  (sql/create-table :files
                    [:analysis_id :text]
                    [:name :text]
                    [:location :text]))

(defn prepare-web-db
  "Prepare input database for storing user and file information in SQLite."
  [db-file]
  (when-not (fs/exists? db-file)
    (sql/with-connection (get-db db-file :create true)
      (create-user-tables)))
  db-file)

(defn get-analyses
  "Retrieve list of analyses run for a specific user."
  [username db-file])

(defn add-analysis
  "Add an analysis and associated files to the database."
  [username analysis db-file]
  (println username)
  (println analysis))