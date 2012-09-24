(ns bcbio.variation.custom.core
  "Commandline dispatch for custom one-off exploratory code"
  (:require [bcbio.variation.custom.nist :as nist]))

(defn -main [prog & args]
  (apply (case (keyword prog)
           :nist nist/summarize-discordants)
         args))