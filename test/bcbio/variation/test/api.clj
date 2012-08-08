(ns bcbio.variation.test.api
  "Tests for top level API for variant metrics"
  (:use [midje.sweet]
        [bcbio.variation.api.metrics])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.variation.index.metrics :as im]))

(background
 (around :facts
         (let [data-dir (str (fs/file "." "test" "data"))
               ref (str (fs/file data-dir "GRCh37.fa"))
               vcf1 (str (fs/file data-dir "gatk-calls.vcf"))
               out-index (str (itx/file-root vcf1) "-metrics.db")]
           (doseq [x [out-index]]
             (itx/remove-path x))
           ?form
           )))

(facts "Retrieve available metrics from a variant file"
  (map :id (available-metrics vcf1 {} nil)) => (has-prefix ["QUAL" "DP"])
  (-> (get-raw-metrics vcf1 ref) first keys) => (contains ["QUAL" :id])
  (let [out (plot-ready-metrics vcf1 ref)]
    (:filename out) => vcf1
    (-> out :metrics first :id) => "QUAL"
    (-> out :metrics first :x-scale :domain) => (just [0.0 10000.0])
    (apply + (-> out :metrics first :vals)) => 1.0))

(facts "Index file for rapid retrieval of variants."
  (im/index-variant-file vcf1 ref) => out-index
  (-> (im/get-raw-metrics vcf1 ref) first keys) => ["PL" "AD" "HaplotypeScore" "QD" "MQ" "DP" "QUAL" :id])