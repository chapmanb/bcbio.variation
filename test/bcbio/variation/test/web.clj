(ns bcbio.variation.test.web
  "Test code supporting a web-based interface for running scoring."
  (:use [clojure.java.io]
        [midje.sweet]
        [bcbio.variation.web.process])
  (:require [fs.core :as fs]
            [clj-yaml.core :as yaml]
            [bcbio.run.itx :as itx]))

(let [conf-file (str (fs/file "config" "web-processing.yaml"))
      config (-> conf-file slurp yaml/parse-string)
      test-dir (str (fs/file "test" "data"))
      work-info {:id "web-test" :dir (str (fs/file test-dir "web-test"))
                 :comparison-genome "NA00001"
                 :in-files {:variant-file (str (fs/file test-dir "phasing-contestant.vcf"))
                            :region-file (str (fs/file test-dir "phasing-contestant-regions.bed"))}}]
  (against-background [(before :facts (vec (map itx/remove-path [(:dir work-info)])))]
    (facts "Setup working directory for processing from configuration info."
      (create-work-config work-info config) => (has-suffix "process.yaml"))))
