(ns bcbio.variation.filter.train
  "Extract training cases from comparisons for machine learning approaches.
  Based on a comparison, identified potential true positives, false positives
  and false negatives to further tweak classifiers."
  (:use [clojure.java.io]
        [bcbio.variation.multiple :only [prep-cmp-name-lookup]])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]))

(defn- get-train-cases
  "Retrieve cases to use for preparing training sets from supplied inputs.
   Prepares list of maps with :target :cmp :c-files, where the latter contains
   all original comparison files."
  [cmps-orig train-info]
  (letfn [(get-train-case [cmps target cmp]
            {:target target :cmp cmp
             :c-files (:c-files (get cmps [target cmp] (get cmps [cmp target])))})]
    (let [cmps (prep-cmp-name-lookup cmps-orig)]
      (map (partial get-train-case cmps (:target train-info)) (:cmps train-info)))))

(defn- prep-false-negatives
  "Retrieve potential false negatives, discordant calls found in all of the
   comparison cases but not in the target."
  [cases out-base]
  (let [out-file (str out-base "-potential-fns.vcf")]
    (when (itx/needs-run? out-file))
    out-file))

(defn- prep-false-positives
  "Retrieve potential false positives, discordant calls from the target not
   found in any of the comparison cases."
  [cases out-base]
  (let [out-file (str out-base "-potential-fps.vcf")]
    (when (itx/needs-run? out-file))
    out-file))

(defn extract-train-cases
  "Prepare exploratory training cases based on specified inputs"
  [cmps train-info exp config]
  (let [out-dir (str (file (get-in config [:dir :prep] (get-in config [:dir :out])) "train"))
        cases (get-train-cases cmps train-info)
        out-base (str (file out-dir (format "%s-%s" (:sample exp) (:target train-info))))]
    (when (not (fs/exists? out-dir))
      (fs/mkdirs out-dir))
    {:fns (prep-false-negatives cases out-base)
     :fps (prep-false-positives cases out-base)}))