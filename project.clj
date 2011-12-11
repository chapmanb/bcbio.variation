(defproject bcbio.variation "0.0.1"
  :description "Manipluate VCF files using the GATK library."
  :dependencies [[org.clojure/clojure "1.3.0"]
                 [org.clojars.chapmanb/gatk "1.3"]
                 [fs "0.11.0"]]
  :dev-dependencies [[midje "1.3-alpha4"]]
  :java-source-path "src/java"
  :aot [bcbio.variation.vcfwalker bcbio.variation.core]
  :main bcbio.variation.core)
