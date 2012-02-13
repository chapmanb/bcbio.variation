(ns bcbio.variation.web.process
  "Run scoring analysis, handling preparation of input files and run configuration."
  (:import [java.util.UUID])
  (:use [clojure.java.io]
        [ring.util.response :only (redirect)])
  (:require [fs.core :as fs]))

(defn run-scoring
  "Prepare output directory for running scoring analysis, with form inputs."
  [request]
  (letfn [(prep-tmp-dir [request]
            (let [tmp-dir (-> request :config :dir :work)
                  cur-dir (fs/file tmp-dir (str (java.util.UUID/randomUUID)))]
              (fs/mkdirs cur-dir)
              (str cur-dir)))
          (download-file [tmp-dir request name]
            (let [cur-param (-> request :params (get name))
                  out-file (fs/file tmp-dir (:filename cur-param))]
              [(keyword name) (if (> (:size cur-param) 0)
                                (do
                                  (copy (:tempfile cur-param) out-file)
                                  (str out-file)))]))]
    (let [tmp-dir (prep-tmp-dir request)
          in-files (into {} (map (partial download-file tmp-dir request)
                                 ["variant-file" "region-file"]))]
      (println tmp-dir)
      (println in-files)))
  (redirect "/"))
