(ns bcbio.variation.validate
  "Combine calls from a multiple technology comparison to produce a set of final
  variants plus a list for validation. The inputs are:
   - Target technology: The outlier technology for picking additional targets. This
     should be well understood enough to set threshold for validation.
   - Validation info: Details for prepping a set of variants for validation
      - thresholds: min and max thresholds for validation
      - approach: validation along the full range of thresholds, or validate top variants
      - count: total number of variants for validation.
  Produces:
   - Final calls
       - calls that overlap in all of the technologies
       - calls that overlap in all but the target, where the target technology quality
         is below the validation threshold.
   - Validate calls
       - calls that overlap in all but the target and fall below configurable threshold.
         These are either sampled from the distribution or picked off the top."
  (:use [ordered.map :only [ordered-map]]
        [bcbio.variation.multiple :only [prep-cmp-name-lookup
                                         multiple-overlap-analysis]]))

(defn get-final-and-tovalidate
  "Prepare files of calls: finalized and validation targets."
  [cmps finalizer config]
  (let [cmps-by-name (prep-cmp-name-lookup (vals cmps) :remove-mods? true
                                           :ignore #{"all" "validate"})
        multi-prep (multiple-overlap-analysis cmps-by-name config (:target finalizer)
                                              :dirname "validate")]
    (ordered-map :final (:true-positives multi-prep)
                 :validate (:false-negatives multi-prep))))

(defn pipeline-validate
  "High level pipeline entry for producing final and to-validate call sets."
  [cmps finalizer exp config]
  {:c-files (get-final-and-tovalidate cmps finalizer config)
   :c1 {:name (:target finalizer)}
   :c2 {:name "validate"}
   :exp exp :dir (:dir config)})
