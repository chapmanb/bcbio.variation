(ns bcbio.variation.core
  (:import [org.broadinstitute.sting.gatk CommandLineGATK])
  (:require [clojure.string :as string]
            [bcbio.align.reorder]
            [bcbio.variation.combine]
            [bcbio.variation.compare]
            [bcbio.variation.ensemble]
            [bcbio.variation.filter.custom]
            [bcbio.variation.haploid]
            [bcbio.variation.utils.core])
  (:gen-class))

(def ^{:doc "Mapping of special command line arguments to main functions"
       :private true}
  altmain-map
  {:compare bcbio.variation.compare/-main
   :ensemble bcbio.variation.ensemble/-main
   :filter bcbio.variation.filter.custom/-main
   :haploid bcbio.variation.haploid/-main
   :prep bcbio.variation.combine/-main
   :reorder bcbio.align.reorder/-main
   :utils bcbio.variation.utils.core/-main})

(defn- get-altmain-fn
  "Retrieve alternative main functions based on first argument."
  [arg]
  (when (and (not (nil? arg))
             (.startsWith arg "variant-"))
    (get altmain-map
         (keyword (string/replace-first arg "variant-" "")))))

(defn -main [& args]
  (if-let [alt-fn (get-altmain-fn (first args))]
    (do
      (apply alt-fn (rest args))
      (shutdown-agents)
      (System/exit 0))
    (CommandLineGATK/main (into-array (if-not (nil? args) args ["-h"])))))
