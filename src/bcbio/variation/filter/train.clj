(ns bcbio.variation.filter.train
  "Extract training cases from comparisons for machine learning approaches.
  Based on a comparison, identified potential true positives, false positives
  and false negatives to further tweak classifiers.")

(defn extract-train-cases
  "Prepare exploratory training cases based on specified inputs"
  [cmps])