(defproject bcbio.variation "0.0.1-SNAPSHOT"
  :description "Clojure API for variation data, built on GATK"
  :dependencies [[org.clojure/clojure "1.3.0"]
                 [org.clojure/math.combinatorics "0.0.2"]
                 [org.clojars.chapmanb/gatk "1.4"]
                 [org.clojars.chapmanb/picard "1.58"]
                 [incanter/incanter-core "1.3.0-SNAPSHOT"]
                 [incanter/incanter-charts "1.3.0-SNAPSHOT"]
                 [fs "1.0.0"]
                 [clj-yaml "0.3.1"]
                 [doric "0.7.0-SNAPSHOT"]
                 [ordered "1.0.0"]]
  :dev-dependencies [[midje "1.3.0"]]
  :java-source-path "src/java"
  :omit-source true
  :aot [bcbio.variation.vcfwalker bcbio.variation.core]
  :main bcbio.variation.core
  :run-aliases {:compare bcbio.variation.compare})
