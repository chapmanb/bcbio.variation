(ns bcbio.variation.core
  (:import [org.broadinstitute.sting.gatk CommandLineGATK])
  (:gen-class))

(defn -main [& args]
  (CommandLineGATK/main (into-array (if-not (nil? args) args ["-h"]))))
