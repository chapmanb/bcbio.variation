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
         These are either sampled from the distribution or picked off the top.")

(defn get-final-and-tovalidate
  "Prepare files of calls: finalized and those that require validation."
  [cmps finalizer config]
  (println finalizer)
  (println (keys cmps))
  (throw (Exception.)))

(defn pipeline-validate
  "High level pipeline entry for producing final and to-validate call sets."
  [cmps finalizer exp config]
  {:c-files (get-final-and-tovalidate cmps finalizer config)
   :c1 {:name (:target finalizer)}
   :c2 {:name "final"}
   :exp exp :dir (:dir config)})
