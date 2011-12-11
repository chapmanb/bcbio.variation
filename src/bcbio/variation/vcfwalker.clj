;; Simple walker to parse a VCF file.

(ns bcbio.variation.vcfwalker
  (:import [bcbio.variation BaseVariantWalker])
  (:use [bcbio.variation.variantcontext :only [from-vc]])
  (:gen-class
   :name bcbio.variation.vcfwalker.VcfSimpleStatsWalker
   :extends bcbio.variation.BaseVariantWalker))

(defn -map [this tracker ref context]
  (if-not (nil? tracker)
    (doseq [vc (map from-vc
                    (.getValues tracker (.variants (.invrns this))
                                (.getLocation context)))]
      (println (:chr vc) (:start vc) (get (:attributes vc) "DP")))))

(defn -reduceInit [this])

(defn -reduce [this value sum])

(defn -onTravelsalDone [this result])
