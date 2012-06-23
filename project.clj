(defproject bcbio.variation "0.0.1-SNAPSHOT"
  :description "Clojure API for variation data, built on GATK"
  :dependencies [[org.clojure/clojure "1.4.0"]
                 [org.clojure/math.combinatorics "0.0.2" :exclusions [org.clojure/clojure]]
                 [org.clojure/data.csv "0.1.2" :exclusions [org.clojure/clojure]]
                 [org.clojars.chapmanb/gatk "1.6.5"]
                 [org.clojars.chapmanb/picard "1.64"]
                 [clj-genomespace "0.1-SNAPSHOT"]
                 [incanter/incanter-core "1.3.0-SNAPSHOT" :exclusions [org.clojure/clojure]]
                 [incanter/incanter-charts "1.3.0-SNAPSHOT" :exclusions [org.clojure/clojure]]
                 [nz.ac.waikato.cms.weka/weka-stable "3.6.6"]
                 [org.clojars.chapmanb/fast-random-forest "0.98"]
                 [com.leadtune/clj-ml "0.2.2" :exclusions [lt/weka hr.irb/fastRandomForest
                                                           org.clojure/clojure]]
                 [fs "1.1.2" :exclusions [org.clojure/clojure]]
                 [clj-yaml "0.3.1"]
                 [doric "0.7.0" :exclusions [org.clojure/clojure]]
                 [ordered "1.0.0" :exclusions [org.clojure/clojure]]
                 [de.kotka/lazymap "3.0.0"]
                 [pallet-fsm "0.1.0"]
                 [clj-logging-config "1.9.7"]
                 [clj-time "0.4.3"]
                 [clj-aws-s3 "0.3.1"]
                 [noir "1.2.2" :exclusions [org.clojure/clojure]]
                 [ring-anti-forgery "0.1.3"]
                 [fetch "0.1.0-alpha2"]
                 [crate "0.2.0-alpha3"]
                 [enlive "1.0.0"]
                 [hiccup "0.3.8"]
                 [domina "1.0.0-beta4" :exclusions [org.clojure/clojurescript]]
                 [jayq "0.1.0-alpha4"]
                 [com.keminglabs/chosen "0.1.6"]]
  :dev-dependencies [[midje "1.4.0" :exclusions [org.clojure/clojure]]
                     [lein-midje "1.0.10"]]
  :plugins [[lein-cljsbuild "0.2.1"]
            [lein-marginalia "0.7.1"]]
  :java-source-path "src/java"
  ;:jvm-opts ["-Xmx4g"]
  :omit-source true
  :aot [bcbio.variation.vcfwalker bcbio.variation.core bcbio.variation.annotate.nbq]
  :main bcbio.variation.core
  :run-aliases {:compare bcbio.variation.compare
                :web bcbio.variation.web.server
                :evaluate bcbio.variation.evaluate
                :gms bcbio.variation.utils.gms
                :haploid bcbio.variation.haploid
                :background bcbio.variation.utils.background
                :popfreq bcbio.variation.utils.popfreq
                :summarize bcbio.variation.utils.summarize
                :reorder bcbio.align.reorder
                :vctest bcbio.variation.variantcontext}
  :cljsbuild {:builds
              [{:builds nil,
                :source-path "src/cljs",
                :compiler
                {:output-to "public/js/score.js",
                 :optimizations :advanced
                 :externs ["externs/jquery.js"]
                 :pretty-print false}}]})
