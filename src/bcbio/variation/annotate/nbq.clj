(ns bcbio.variation.annotate.nbq
  "GATK annotator that calculates Mean Neighboring Base Quality (NBQ) for variants.

  The motivation for this annotation is that regional base quality influences whether
  a call is correct. The Atlas2 paper describes the metric in more detail:

  http://www.biomedcentral.com/1471-2105/13/8/abstract
  "
  (:import [org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation]
           [org.broadinstitute.variant.vcf VCFInfoHeaderLine VCFHeaderLineType]
           [org.broadinstitute.sting.utils BaseUtils])
  (:require [criterium.stats :as stats])
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
    - Filter reads to those that match an alternative base
    - Get quality from reads and pull out qualities in surrounding region
    - Calculate mean and return."
  [_ _ _ _ contexts vc _]
  (letfn [(orient-reads [[offset read]]
            (if (.getReadNegativeStrandFlag read)
              {:offset offset
               :bases (BaseUtils/simpleReverseComplement (.getReadBases read))
               :quals (-> read .getBaseQualities vec reverse)}
              {:offset offset
               :bases (.getReadString read)
               :quals (-> read .getBaseQualities vec)}))
          (neighbor-qualities [{:keys [offset quals]}]
            (map #(nth quals % nil) (range (- offset flank-bp) (+ offset flank-bp))))
          (supports-alt? [alt-bases {:keys [offset bases]}]
            (let [base (char (nth bases offset))]
              (contains? alt-bases base)))
          (pileup-qualities [alt-bases pileup]
            (->> (map vector (.getOffsets pileup) (.getReads pileup))
                 (map orient-reads)
                 (filter (partial supports-alt? alt-bases))
                 (map neighbor-qualities)))]
    (let [alt-bases (->> (.getAlternateAlleles vc)
                              (map #(.getBaseString %))
                              (map first)
                              set)]
      (when (seq alt-bases)
        {"NBQ" (->> contexts
                    vals
                    (map #(.getBasePileup %))
                    (map (partial pileup-qualities alt-bases))
                    flatten
                    (remove nil?)
                    stats/mean
                    float
                    (format "%.2f"))}))))
