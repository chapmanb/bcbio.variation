(ns bcbio.variation.test.web
  "Test code supporting a web-based interface for running scoring."
  (:use [clojure.java.io]
        [midje.sweet]
        [bcbio.variation.web.process])
  (:require [fs.core :as fs]
            [clj-yaml.core :as yaml]))

(defn remove-path [x]
  (if (fs/exists? x)
    (if (fs/directory? x)
      (fs/delete-dir x)
      (fs/delete x))))

(let [conf-file (str (fs/file "config" "web-processing.yaml"))
      config (-> conf-file slurp yaml/parse-string)
      test-dir (str (fs/file "test" "data"))
      work-info {:id "web-test" :dir (str (fs/file test-dir "web-test"))}
      in-files {:variant-file (str (fs/file test-dir "phasing-contestant.vcf"))
                :region-file (str (fs/file test-dir "phasing-contestant-regions.bed"))}]
  (against-background [(before :facts (vec (map remove-path [(:dir work-info)])))]
    (facts "Setup working directory for processing from configuration info."
      (create-work-config in-files work-info config) => (has-suffix "process.yaml"))))
