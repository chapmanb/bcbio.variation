(defproject bcbio.variation "0.2.2"
  :description "Toolkit to analyze genomic variation data, built on the GATK with Clojure"
  :license {:name "MIT" :url "http://www.opensource.org/licenses/mit-license.html"}
  :dependencies [[org.clojure/clojure "1.5.1"]
                 [org.clojure/math.combinatorics "0.0.3" :exclusions [org.clojure/clojure]]
                 [org.clojure/data.csv "0.1.2" :exclusions [org.clojure/clojure]]
                 [org.clojure/tools.cli "0.2.2"]
                 [clj-stacktrace "0.2.5"]
                 [bcbio.run "0.0.1"]
                 ;; GATK requirements
                 [org.clojars.chapmanb/gatk-engine "3.2"]
                 [org.clojars.chapmanb/gatk-tools-public "3.2"]
                 [org.clojars.chapmanb/gatk-utils "3.2"]
                 [org.clojars.chapmanb/picard "1.112"]
                 [org.clojars.chapmanb/htsjdk "1.112"]
                 [colt/colt "1.2.0"]
                 [commons-lang "2.5"]
                 [log4j "1.2.17"]
                 [org.slf4j/slf4j-log4j12 "1.7.7"]
                 [org.apache.commons/commons-jexl "2.1.1"]
                 [org.reflections/reflections "0.9.9-RC1"]
                 [org.simpleframework/simple-xml "2.0.4"]
                 [org.apache.servicemix.bundles/org.apache.servicemix.bundles.jets3t "0.8.1_1"]
                 ;;
                 [org.biojava/biojava3-core "4.0.0-SNAPSHOT"]
                 [org.biojava/biojava3-alignment "4.0.0-SNAPSHOT"]
                 [org.clojars.chapmanb/circdesigna "0.0.2" :exclusions [net.sf.beaver/beaver-ant]]
                 [nz.ac.waikato.cms.weka/weka-stable "3.6.6"]
                 [org.clojars.chapmanb/fast-random-forest "0.98"]
                 [com.leadtune/clj-ml "0.2.4" :exclusions [cs.waikato.ac.nz/weka
                                                           hr.irb/fastRandomForest
                                                           org.clojure/clojure
                                                           incanter/incanter-core
                                                           incanter/incanter-charts]]
                 [clj-yaml "0.4.0"]
                 [doric "0.8.0" :exclusions [org.clojure/clojure]]
                 [iota "1.1.1"]
                 [ordered "1.3.2" :exclusions [org.clojure/clojure]]
                 [de.kotka/lazymap "3.1.1"]
                 [lonocloud/synthread "1.0.4"]
                 [pallet-fsm "0.1.0"]
                 [criterium "0.4.1" :exclusions [org.clojure/clojure]]
                 [clj-time "0.5.1"]
                 [org.clojure/java.jdbc "0.2.2"]
                 [org.xerial/sqlite-jdbc "3.7.2"]
                 [c3p0/c3p0 "0.9.1.2"]
                 [hiccup "1.0.1"]
                 [enlive "1.0.1" :exclusions [org.clojure/clojure]]]
  :min-lein-version "2.0.0"
  :plugins [[lein-marginalia "0.7.1"]
            [lein-midje "3.0.1"]]
  :profiles {:dev {:dependencies
                   ;; Testing dependencies
                   [[midje "1.5.1" :exclusions [org.clojure/clojure ordered]]
                    ;; Non-uberjar dependencies we should extract into separate functionality
                    [incanter/incanter-core "1.5.1" :exclusions [org.clojure/clojure junit]]
                    [incanter/incanter-charts "1.5.1" :exclusions [org.clojure/clojure junit]]
                    [incanter/incanter-excel "1.5.1" :exclusions [org.clojure/clojure junit]]
                    [org.clojars.chapmanb/snpeff "3.1" :exclusions [com.googlecode.charts4j/charts4j
                                                                    junit]]
                    [clj-genomespace "0.1.3"]
                    [clj-blend "0.1.1-SNAPSHOT"]
                    [clj-aws-s3 "0.3.1" :exclusions [org.codehaus.jackson/jackson-mapper-asl
                                                     org.codehaus.jackson/jackson-core-asl]]]}}
  :repositories {"sonatype-snapshots" {:url "http://oss.sonatype.org/content/repositories/snapshots"
                                       :snapshots true}}
  :java-source-paths ["src/java"]
  :javac-options ["-nowarn" "-target" "1.6" "-source" "1.6"]
  ;;:jvm-opts ["-Xms750m" "-Xmx2g"]
  :omit-source false
  :aot [bcbio.variation.vcfwalker bcbio.variation.core bcbio.variation.annotate.nbq
        bcbio.variation.annotate.entropy bcbio.variation.annotate.mfe]
  :main bcbio.variation.core
  :aliases {"variant-compare" ["run" "-m" "bcbio.variation.compare"]
            "variant-prep" ["run" "-m" "bcbio.variation.combine"]
            "variant-evaluate" ["run" "-m" "bcbio.variation.evaluate"]
            "variant-ensemble" ["run" "-m" "bcbio.variation.ensemble"]
            "variant-haploid" ["run" "-m" "bcbio.variation.haploid"]
            "variant-recall" ["run" "-m" "bcbio.variation.recall"]
            "variant-reorder" ["run" "-m" "bcbio.align.reorder"]
            "variant-utils" ["run" "-m" "bcbio.variation.utils.core"]
            "variant-custom" ["run" "-m" "bcbio.variation.custom.core"]})
