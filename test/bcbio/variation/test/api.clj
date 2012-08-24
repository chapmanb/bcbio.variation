(ns bcbio.variation.test.api
  "Tests for top level API for variant metrics"
  (:use [midje.sweet]
        [bcbio.variation.api.metrics]
        [bcbio.variation.api.shared :only [set-config-from-file!]])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.variation.index.gemini :as gemini]
            [bcbio.variation.index.metrics :as im]))

(background
 (around :facts
         (let [data-dir (str (fs/file "." "test" "data"))
               web-yaml (str (fs/file "." "config" "web-processing.yaml"))
               ref (str (fs/file data-dir "GRCh37.fa"))
               vcf1 (str (fs/file data-dir "gatk-calls.vcf"))
               out-index (str (itx/file-root vcf1) "-metrics.db")
               gemini-index (str (itx/file-root vcf1) "-gemini.db")]
           (doseq [x [out-index gemini-index]]
             (itx/remove-path x))
           ?form
           )))

(facts "Retrieve available metrics from a variant file"
  (set-config-from-file! web-yaml)
  (map :id (available-metrics vcf1)) => (has-prefix ["QUAL" "DP"])
  (-> (get-raw-metrics vcf1) first keys) => (contains ["QUAL" "DP" :id] :in-any-order :gaps-ok)
  (let [out (plot-ready-metrics vcf1)]
    (:filename out) => vcf1
    (-> out :metrics first :id) => "QUAL"
    (-> out :metrics first :x-scale :domain) => (just [0.0 10000.0])
    (apply + (-> out :metrics first :vals)) => (roughly 1.0)))

(facts "Index file for rapid retrieval of variants."
  (im/index-variant-file vcf1 ref) => out-index
  (-> (im/get-raw-metrics vcf1 ref) first keys) => ["PL" "AD" "HaplotypeScore" "QD" "MQ" "DP" "QUAL" :id])

(facts "Index and retrieve metrics using Gemini."
  (when-let [idx (gemini/index-variant-file vcf1 ref)]
    idx => gemini-index
    (-> (gemini/get-raw-metrics vcf1 ref) first keys) => (contains [:id "gms_illumina" "aaf_1kg_all"]
                                                                   :in-any-order)))