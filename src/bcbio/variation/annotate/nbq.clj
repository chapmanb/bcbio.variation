(ns bcbio.variation.annotate.nbq
  "GATK annotator that calculates Mean Neighboring Base Quality (NBQ) for variants.

  The motivation for this annotation is that regional base quality influences whether
  a call is correct. The Atlas2 paper describes the metric in more detail:

  http://www.biomedcentral.com/1471-2105/13/8/abstract
  "
  (:import [org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation]
           [org.broadinstitute.sting.utils.codecs.vcf VCFInfoHeaderLine VCFHeaderLineType])
  (:require [incanter.stats :as istats])
  (:gen-class
   :name bcbio.variation.annotate.nbq.MeanNeighboringBaseQuality
   :extends org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation))

(def flank-bp 5)

(defn -getKeyNames
  [_]
  ["NBQ"])

(defn -getDescriptions
  [_]
  [(VCFInfoHeaderLine. "NBQ" 1 VCFHeaderLineType/Float
                       (format "Mean Neighboring Base Quality, includes %sbp on both sides"
                               flank-bp))])

(defn -annotate
  "Provide Mean Neighboring Base Quality calculations at a position.

    - Get a pileup for each sample context.
    - Use pileup to retrieve reads and current offsets.
    - Get quality from reads and pull out qualities in surrounding region
    - Calculate mean and return."
  [_ _ _ _ contexts _]
  (letfn [(get-pileup [context]
            (if (.hasExtendedEventPileup context)
              (.getExtendedEventPileup context)
              (.getBasePileup context)))
          (neighbor-qualities [[offset read]]
            (let [quals (-> read .getBaseQualities vec)]
              (map #(nth quals % nil) (range (- offset flank-bp) (+ offset flank-bp)))))
          (pileup-qualities [pileup]
            (map neighbor-qualities (map vector (.getOffsets pileup) (.getReads pileup))))]
    {"NBQ" (->> contexts
                vals
                (map get-pileup)
                (map pileup-qualities)
                flatten
                (remove nil?)
                istats/mean
                (format "%.2f"))}))
