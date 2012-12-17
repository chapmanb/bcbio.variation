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
        [bcbio.variation.combine :only [combine-variants]]
        [bcbio.variation.multiple :only [prep-cmp-name-lookup
                                         multiple-overlap-analysis]]
        [bcbio.variation.report :only [count-variants]]
        [bcbio.variation.variantcontext :only [parse-vcf get-vcf-iterator
                                               write-vcf-w-template]])
  (:require [bcbio.run.broad :as broad]
            [bcbio.run.itx :as itx]))

(defn- select-by-general
  "Base functionality for subsetting a file with SelectVariants."
  [select-args ext in-vcf ref]
  (let [file-info {:out-vcf (itx/add-file-part in-vcf ext)}
        args (concat ["-R" ref
                      "--variant" in-vcf
                      "-o" :out-vcf]
                      select-args)]
    (broad/run-gatk "SelectVariants" args file-info {:out [:out-vcf]})
    (:out-vcf file-info)))

(defn select-by-filters
  "Subset a VCF file with specific hard filters."
  [filters in-vcf ext ref]
  (select-by-general (interleave (repeat "--select_expressions") filters)
                     ext in-vcf ref))

;; ## Validation targets by random sampling

(defn select-by-random
  "Subset a VCF file with a random number of variants."
  [n in-vcf ref]
  (let [total (count-variants in-vcf ref (fn [x] true)) 
        frac (if (pos? total) (float (/ n total)) 0.0)]
    (select-by-general ["--select_random_fraction" frac] "randsubset" in-vcf ref)))

(defmulti get-to-validate
  "Select set of variants to validate from total set of potentials."
  (fn [in-vcf finalizer ref] (keyword (get-in finalizer [:params :validate :approach]))))

(defmethod get-to-validate :random
  [in-vcf finalizer ref]
  (select-by-random (get-in finalizer [:params :validate :count]) in-vcf ref))

; ## Validation targets by sorting

(defn- extract-sort-metrics
  "Provide function to extract metric used in sorting from a variant context.
  Returns a list with the first being the count of items found in set overlap
  and the second the metric to sort by. Currently only handles a single metric
  for sorting."
  [finalizer]
  {:pre [(= 1 (count (get-in finalizer [:params :validate :top-metric])))]}
  (let [metric (get-in finalizer [:params :validate :top-metric 0 :name])
        mod (get-in finalizer [:params :validate :top-metric 0 :mod])]
    (fn [vc]
      (let [base (-> vc :attributes (get metric "-1000.0") (Float/parseFloat))]
        [(count (re-seq (re-pattern (:target finalizer))
                           (-> vc :attributes (get "set" ""))))
         (if mod (* mod base) base)]))))

(defn- get-top-variants
  "Retrieve top variants sorted by metrics of interest."
  [vcf-file finalizer ref]
  (with-open [vcf-iter (get-vcf-iterator vcf-file ref)]
    (let [metric-gettr (extract-sort-metrics finalizer)]
      (set
       (map (juxt :chr :start)
            (take (get-in finalizer [:params :validate :count])
                  (reverse
                   (sort-by metric-gettr (parse-vcf vcf-iter)))))))))

(defmethod get-to-validate :top
  [in-vcf finalizer ref]
  (let [out-file (itx/add-file-part in-vcf "topsubset")]
    (when (itx/needs-run? out-file)
      (let [to-keep (get-top-variants in-vcf finalizer ref)]
        (with-open [vcf-iter (get-vcf-iterator in-vcf ref)]
          (write-vcf-w-template in-vcf {:out out-file}
                                (map :vc
                                     (filter #(contains? to-keep ((juxt :chr :start) %))
                                             (parse-vcf vcf-iter)))
                                ref))))
    out-file))

(defn get-final-and-tovalidate
  "Prepare files of calls: finalized and validation targets."
  [cmps finalizer config]
  (let [cmps-by-name (prep-cmp-name-lookup (vals cmps) :remove-mods? true
                                           :ignore #{"all" "validate"})
        ref (-> cmps-by-name vals first :exp :ref)
        multi-prep (multiple-overlap-analysis cmps-by-name config (:target finalizer)
                                              :dirname "validate")]
    (ordered-map
     :final (if-let [keep-filters (get-in finalizer [:params :filters :keep])]
              (combine-variants [(:true-positives multi-prep)
                                 (select-by-filters keep-filters (:false-negatives multi-prep)
                                                    "keepsubset" ref)]
                                ref :merge-type :full)
              (:true-positives multi-prep))
     :validate (get-to-validate
                (let [orig (:target-overlaps multi-prep)]
                  (if-let [val-filters (get-in finalizer [:params :filters :validate])]
                    (select-by-filters val-filters orig "checksubset" ref)
                    orig))
                finalizer ref))))

(defn pipeline-validate
  "High level pipeline entry for producing final and to-validate call sets."
  [cmps finalizer exp config]
  {:c-files (get-final-and-tovalidate cmps finalizer config)
   :c1 {:name (:target finalizer)}
   :c2 {:name "validate"}
   :exp exp :dir (:dir config)})
