(defproject bcbio.variation "0.0.1-SNAPSHOT"
  :description "Clojure API for variation data, built on GATK"
  :dependencies [[org.clojure/clojure "1.4.0"]
                 [org.clojure/math.combinatorics "0.0.2" :exclusions [org.clojure/clojure]]
                 [org.clojure/data.csv "0.1.2" :exclusions [org.clojure/clojure]]
                 ;; GATK requirements
                 [org.clojars.chapmanb/gatk-lite "2.0.21"]
                 [org.clojars.chapmanb/picard "1.73"]
                 [org.clojars.chapmanb/sam "1.73"]
                 [org.clojars.chapmanb/tribble "110"]
                 [org.clojars.chapmanb/cofoja "1.0-20110609"]
                 [org.clojars.chapmanb/jama "1.0.2"]
                 [org.apache.commons/commons-jexl "2.1.1"]
                 [org.apache.commons/commons-math "2.2"]
                 [org.reflections/reflections "0.9.5-RC2"]
                 [org.simpleframework/simple-xml "2.0.4"]
                 [colt/colt "1.2.0"]
                 ;;
                 [org.biojava/biojava3-core "3.0.4"]
                 [org.biojava/biojava3-alignment "3.0.4"]
                 [clj-genomespace "0.1.2-SNAPSHOT"]
                 [incanter/incanter-core "1.3.0" :exclusions [org.clojure/clojure]]
                 [incanter/incanter-charts "1.3.0" :exclusions [org.clojure/clojure]]
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
                 [clj-time "0.4.3"]
                 [clj-aws-s3 "0.3.1"]
                 [org.clojure/java.jdbc "0.2.2"]
                 [org.xerial/sqlite-jdbc "3.7.2"]
                 [noir "1.2.2" :exclusions [org.clojure/clojure]]
                 [ring-anti-forgery "0.1.3"]
                 [fetch "0.1.0-alpha2" :exclusions [org.clojure/clojure]]
                 [crate "0.2.0-alpha3"]
                 [enlive "1.0.0"]
                 [hiccup "0.3.8"]
                 [domina "1.0.0-beta4" :exclusions [org.clojure/clojurescript]]
                 [jayq "0.1.0-alpha4"]
                 [com.keminglabs/chosen "0.1.6"]]
  :profiles {:dev
             {:dependencies
              [[midje "1.4.0" :exclusions [org.clojure/clojure]]]}}
  :plugins [[lein-cljsbuild "0.2.1"]
            [lein-marginalia "0.7.1"]
            [lein-midje "2.0.0-SNAPSHOT"]]
  :repositories {"biojava" {:url "http://www.biojava.org/download/maven/"
                            :snapshots false}}
  :java-source-paths ["src/java"]
 ;:jvm-opts ["-Xmx4g"]
  :omit-source true
  :aot [bcbio.variation.vcfwalker bcbio.variation.core bcbio.variation.annotate.nbq]
  :main bcbio.variation.core
  :aliases {"variant-compare" ["run" "-m" "bcbio.variation.compare"]
            "variant-web" ["run" "-m" "bcbio.variation.web.server"]
            "variant-evaluate" ["run" "-m" "bcbio.variation.evaluate"]
            "variant-gms" ["run" "-m" "bcbio.variation.utils.gms"]
            "variant-haploid" ["run" "-m" "bcbio.variation.haploid"]
            "variant-utils" ["run" "-m" "bcbio.variation.utils.core"]
            "variant-recall" ["run" "-m" "bcbio.variation.recall"]
            "variant-reorder" ["run" "-m" "bcbio.align.reorder"]
            "variant-vctest" ["run" "-m" "bcbio.variation.variantcontext"]}
  :cljsbuild {:builds
              [{:builds nil
                :source-path "src/cljs"
                :compiler
                {:output-to "public/js/score.js"
                 :optimizations :advanced
                 :externs ["externs/jquery.js"]
                 :pretty-print false}}]})
