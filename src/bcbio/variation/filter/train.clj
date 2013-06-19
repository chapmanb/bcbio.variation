(ns bcbio.variation.filter.train
  "Extract training cases from comparisons for machine learning approaches.
  Based on a comparison, identified potential true positives, false positives
  and false negatives to further tweak classifiers."
  (:use [clojure.java.io]
        [bcbio.variation.multiple :only [prep-cmp-name-lookup]])
  (:require [me.raynes.fs :as fs]
            [bcbio.run.broad :as broad]
            [bcbio.run.itx :as itx]))

(defn- select-concordant
  "Retrieve output file of concordant calls between two sets of variant calls"
  [fname1 fname2 ref out-file]
  (let [args ["-R" ref
              "--variant" fname1
              "--concordance" fname2
              "--out" :out-vcf]]
    (broad/run-gatk "SelectVariants" args {:out-vcf out-file} {:out [:out-vcf]}))
  out-file)

(defn- prep-common
  "Common infrastructure for generating training values"
  [case-kw file-ext cases out-base ref-file]
  (letfn [(get-discordant-by-kw [x]
            (get-in x [:c-files (keyword (str (get x case-kw) "-discordant"))]))]
    (let [out-file (str out-base file-ext)]
      (when (itx/needs-run? out-file)
        (apply select-concordant
               (concat (map get-discordant-by-kw cases) [ref-file out-file])))
      out-file)))

(defn- prep-false-negatives
  "Retrieve potential false negatives, discordant calls found in all of the
   comparison cases but not in the target."
  [cases out-base ref-file]
  (prep-common :cmp "-potential-fns.vcf"
               cases out-base ref-file))

(defn- prep-false-positives
  "Retrieve potential false positives, discordant calls from the target not
   found in any of the comparison cases."
  [cases out-base ref-file]
  (prep-common :target "-potential-fps.vcf"
               cases out-base ref-file))

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

(defn extract-train-cases
  "Prepare exploratory training cases based on specified inputs"
  [cmps train-info exp config]
  (let [out-dir (str (file (get-in config [:dir :prep] (get-in config [:dir :out])) "train"))
        cases (get-train-cases cmps train-info)
        out-base (str (file out-dir (format "%s-%s" (:sample exp) (:target train-info))))]
    (when (not (fs/exists? out-dir))
      (fs/mkdirs out-dir))
    {:fns (prep-false-negatives cases out-base (:ref exp))
     :fps (prep-false-positives cases out-base (:ref exp))}))