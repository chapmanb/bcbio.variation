(ns bcbio.variation.test.api
  "Tests for top level API for variant metrics"
  (:use [midje.sweet]
        [bcbio.variation.api.metrics])
  (:require [fs.core :as fs]))

(background
 (around :facts
         (let [data-dir (str (fs/file "." "test" "data"))
               ref (str (fs/file data-dir "GRCh37.fa"))
               vcf1 (str (fs/file data-dir "gatk-calls.vcf"))]
           ?form
           )))

(facts "Retrieve available metrics from a variant file"
  (map :id (available-metrics vcf1)) => (has-prefix ["QUAL" "DP"])
  (let [out (plot-ready-metrics vcf1 ref)]
    (:filename out) => vcf1
    (-> out :metrics first :id) => "QUAL"
    (-> out :metrics first :x-scale :domain) => (just [(roughly 307.987) 8830.513])
    (apply + (-> out :metrics first :vals)) => 1.0))