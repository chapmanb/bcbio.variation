(defproject bcbio.variation "0.0.1-SNAPSHOT"
  :description "Clojure API for variation data, built on GATK"
  :dependencies [[org.clojure/clojure "1.3.0"]
                 [org.clojure/math.combinatorics "0.0.2"]
                 [org.clojure/data.csv "0.1.2"]
                 [org.clojars.chapmanb/gatk "1.6.5"]
                 [org.clojars.chapmanb/picard "1.64"]
                 [incanter/incanter-core "1.3.0-SNAPSHOT"]
                 [incanter/incanter-charts "1.3.0-SNAPSHOT"]
                 [nz.ac.waikato.cms.weka/weka-stable "3.6.6"]
                 [org.clojars.chapmanb/fast-random-forest "0.98"]
                 [com.leadtune/clj-ml "0.2.2" :exclusions [lt/weka hr.irb/fastRandomForest]]
                 [fs "1.1.2"]
                 [clj-yaml "0.3.1"]
                 [doric "0.7.0-SNAPSHOT"]
                 [ordered "1.0.0"]
                 [clj-aws-s3 "0.3.1"]
                 [compojure "1.0.1"]
                 [ring "1.0.2"]
                 [enlive "1.0.0"]
                 [domina "1.0.0-beta4" :exclusions [org.clojure/clojurescript]]]
  :dev-dependencies [[midje "1.3.0" :exclusions [org.clojure/clojure]]
                     [lein-midje "1.0.10"]]
  :plugins [[lein-cljsbuild "0.1.9"]]
  :java-source-path "src/java"
  :jvm-opts ["-Xmx4g"]
  :omit-source true
  :aot [bcbio.variation.vcfwalker bcbio.variation.core bcbio.variation.annotate.nbq]
  :main bcbio.variation.core
  :run-aliases {:compare bcbio.variation.compare
                :web bcbio.variation.web.server
                :popfreq bcbio.variation.utils.popfreq
                :background bcbio.variation.utils.background
                :reorder bcbio.align.reorder}
  :cljsbuild {:builds
              [{:builds nil,
                :source-path "src/cljs",
                :compiler
                {:output-to "public/js/score.js",
                 :optimizations :advanced,
                 :pretty-print false}}]})
