(ns bcbio.variation.annotate.mfe
  "Calculate delta G Minimum Free Energy for sequence secondary structures.
   Extracts regions surrounding variants and identifies the free energy
   of the most problematic secondary structures. Larger negative free energy
   values are especially stable and problematic."
  (:import [org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation]
           [org.broadinstitute.sting.utils.codecs.vcf VCFInfoHeaderLine VCFHeaderLineType])
  (:use [circdesigna.core :only [min-free-energy]]
        [bcbio.variation.annotate.entropy :only [get-flank-seq]])
  (:gen-class
   :name bcbio.variation.annotate.entropy.MinFreeEnergy
   :extends org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation))

(def flank-bp 15)

(defn -getKeyNames [_]
  ["MFE"])

(defn -getDescriptions [_]
  [(VCFInfoHeaderLine. "MFE" 1 VCFHeaderLineType/Float
                       (format (str "delta G minimum free energy of the most problematic "
                                    "secondary structure +/- %sbp around variant")
                               flank-bp))])

(defn -annotate
  "Retrieve flanking region surrounding variant and calculate MFE."
  [_ _ _ ref _ _]
  {"MFE" (->> (get-flank-seq ref flank-bp)
              min-free-energy
              (format "%.2f"))})