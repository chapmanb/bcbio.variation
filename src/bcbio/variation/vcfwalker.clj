(ns bcbio.variation.vcfwalker
  "Simple walker to parse a VCF file and display distribution of call
  quality scores"
  (:import [bcbio.variation BaseVariantWalker])
  (:use [bcbio.variation.variantcontext :only [from-vc]])
  (:require ;[incanter.charts :as icharts]
            [incanter.core :as icore])
  (:gen-class
   :name bcbio.variation.vcfwalker.VcfSimpleStatsWalker
   :extends bcbio.variation.BaseVariantWalker))

(defn -map
  "Retrieve VariantContexts and extract the variant quality score."
  [this tracker ref context]
  (if-not (nil? tracker)
    (for [vc (map from-vc
                    (.getValues tracker (.variants (.invrns this))
                                (.getLocation context)))]
      (:qual vc))))

(defn -reduceInit
  "Initialize an empty list to collect our quality information"
  [this]
  [])

(defn -reduce
  "Add current quality information to the collected list."
  [this cur coll]
  (if-not (nil? cur)
    (vec (flatten [coll cur]))
    coll))

(defn -onTraversalDone
  "Plot histogram of quality scores."
  [this result]
  (println result))
;  (doto (icharts/histogram result
;                           :x-label "Variant quality"
;                           :nbins 50)
;    (icore/save (.out this) :width 500 :height 400)))
