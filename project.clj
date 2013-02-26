(defproject bcbio.variation "0.0.7-SNAPSHOT"
  :description "Toolkit to analyze genomic variation data, built on the GATK with Clojure"
  :license {:name "MIT" :url "http://www.opensource.org/licenses/mit-license.html"}
  :dependencies [[org.clojure/clojure "1.4.0"]
                 [org.clojure/math.combinatorics "0.0.2" :exclusions [org.clojure/clojure]]
                 [org.clojure/data.csv "0.1.2" :exclusions [org.clojure/clojure]]
                 [org.clojure/core.match "0.2.0-alpha9"]
                 [org.clojure/tools.cli "0.2.2"]
                 [clj-stacktrace "0.2.5"]
                 ;; GATK requirements
                 [org.clojars.chapmanb/gatk-lite "2.3.4"]
                 [org.clojars.chapmanb/picard "1.73"]
                 [org.clojars.chapmanb/sam "1.73"]
                 [org.clojars.chapmanb/tribble "119"]
                 ;; [org.clojars.chapmanb/cofoja "1.0-r139"]
                 [org.clojars.chapmanb/jama "1.0.2"]
                 [org.apache.commons/commons-jexl "2.1.1"]
                 [org.apache.commons/commons-math "2.2"]
                 [org.reflections/reflections "0.9.8"]
                 [org.simpleframework/simple-xml "2.0.4"]
                 [colt/colt "1.2.0"]
                 ;;
                 [org.clojars.chapmanb/snpeff "3.1"]
                 [org.biojava/biojava3-core "3.0.4"]
                 [org.biojava/biojava3-alignment "3.0.4"]
                 [org.clojars.chapmanb/circdesigna "0.0.2"]
                 [clj-genomespace "0.1.3"]
                 [clj-blend "0.1.1-SNAPSHOT"]
                 [incanter/incanter-core "1.4.0" :exclusions [org.clojure/clojure]]
                 [incanter/incanter-charts "1.4.0" :exclusions [org.clojure/clojure]]
                 [incanter/incanter-excel "1.4.0" :exclusions [org.clojure/clojure]]
                 [nz.ac.waikato.cms.weka/weka-stable "3.6.6"]
                 [org.clojars.chapmanb/fast-random-forest "0.98"]
                 [com.leadtune/clj-ml "0.2.4" :exclusions [cs.waikato.ac.nz/weka
                                                           hr.irb/fastRandomForest
                                                           org.clojure/clojure
                                                           incanter/incanter-core
                                                           incanter/incanter-charts]]
                 [fs "1.1.2" :exclusions [org.clojure/clojure]]
                 [clj-yaml "0.3.1"]
                 [doric "0.8.0" :exclusions [org.clojure/clojure]]
                 [ordered "1.3.2" :exclusions [org.clojure/clojure]]
                 [de.kotka/lazymap "3.1.0"]
                 [pallet-fsm "0.1.0"]
                 [clj-time "0.4.3"]
                 [clj-aws-s3 "0.3.1" :exclusions [org.codehaus.jackson/jackson-mapper-asl
                                                  org.codehaus.jackson/jackson-core-asl]]
                 [org.clojure/java.jdbc "0.2.2"]
                 [org.xerial/sqlite-jdbc "3.7.2"]
                 [c3p0/c3p0 "0.9.1.2"]
                 [hiccup "1.0.1"]
                 [enlive "1.0.1" :exclusions [org.clojure/clojure]]]
  :min-lein-version "2.0.0"
  :plugins [[lein-marginalia "0.7.1"]
            [lein-midje "2.0.0-SNAPSHOT"]]
  :profiles {:dev {:dependencies
                   [[midje "1.4.0" :exclusions [org.clojure/clojure ordered]]]}}
  :repositories {"biojava" {:url "http://www.biojava.org/download/maven/"
                            :snapshots false}}
  :java-source-paths ["src/java"]
  :javac-options ["-nowarn" "-target" "1.6" "-source" "1.6"]
  ;:jvm-opts ["-Xmx4g"]
  :omit-source false
  :aot [bcbio.variation.vcfwalker bcbio.variation.core bcbio.variation.annotate.nbq
        bcbio.variation.annotate.entropy bcbio.variation.annotate.mfe]
  :main bcbio.variation.core
  :aliases {"variant-compare" ["run" "-m" "bcbio.variation.compare"]
            "variant-prep" ["run" "-m" "bcbio.variation.combine"]
            "variant-evaluate" ["run" "-m" "bcbio.variation.evaluate"]
            "variant-haploid" ["run" "-m" "bcbio.variation.haploid"]
            "variant-recall" ["run" "-m" "bcbio.variation.recall"]
            "variant-reorder" ["run" "-m" "bcbio.align.reorder"]
            "variant-utils" ["run" "-m" "bcbio.variation.utils.core"]
            "variant-custom" ["run" "-m" "bcbio.variation.custom.core"]})
