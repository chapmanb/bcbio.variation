;; Simple walker to parse a VCF file and display distribution of call
;; quality scores

(ns bcbio.variation.vcfwalker
  (:import [bcbio.variation BaseVariantWalker])
  (:use [bcbio.variation.variantcontext :only [from-vc]])
  (:require [incanter.charts :as icharts]
            [incanter.core :as icore])
  (:gen-class
   :name bcbio.variation.vcfwalker.VcfSimpleStatsWalker
   :extends bcbio.variation.BaseVariantWalker))

(defn -map [this tracker ref context]
  "Retrieve VariantContexts and extract the variant quality score."
  (if-not (nil? tracker)
    (for [vc (map from-vc
                    (.getValues tracker (.variants (.invrns this))
                                (.getLocation context)))]
      (-> vc :genotypes first :qual))))

(defn -reduceInit [this]
  "Initialize an empty list to collect our quality information"
  [])

(defn -reduce [this cur coll]
  "Add current quality information to the collected list."
  (if-not (nil? cur)
    (flatten [coll cur])
    coll))

(defn -onTraversalDone [this result]
  "Plot histogram of quality scores."
  (doto (icharts/histogram result
                           :x-label "Variant quality"
                           :nbins 50)
    (icore/save (.out this) :width 500 :height 400)))
