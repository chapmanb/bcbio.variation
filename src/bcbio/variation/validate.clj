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
        [bcbio.variation.variantcontext :only [parse-vcf get-vcf-source
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

(defn select-by-random
  "Subset a VCF file with a random number of variants."
  [count in-vcf ref]
  (select-by-general ["--select_random_number" count] "randsubset" in-vcf ref))

(defmulti get-to-validate
  "Select set of variants to validate from total set of potentials."
  (fn [in-vcf params ref] (keyword (get-in params [:validate :approach]))))

(defmethod get-to-validate :random
  [in-vcf params ref]
  (select-by-random (get-in params [:validate :count]) in-vcf ref))

(defn- get-top-variants
  "Retrieve top variants sorted by metrics of interest."
  [vcf-iter params])

(defmethod get-to-validate :top
  [in-vcf params ref]
  (let [out-file (itx/add-file-part in-vcf "topsubset")]
    (when (itx/needs-run? out-file)
      (with-open [vcf-source (get-vcf-source in-vcf ref)]
        (write-vcf-w-template in-vcf out-file
                              (get-top-variants (parse-vcf vcf-source) params)
                              ref)))
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
     :final (combine-variants [(:true-positives multi-prep)
                               (select-by-filters (get-in finalizer [:params :filters :keep])
                                                  (:false-negatives multi-prep) "keepsubset" ref)]
                              ref :merge-type :full)
     :validate (get-to-validate (select-by-filters (get-in finalizer [:params :filters :validate])
                                                   (:false-negatives multi-prep) "checksubset" ref)
                                (:params finalizer) ref))))

(defn pipeline-validate
  "High level pipeline entry for producing final and to-validate call sets."
  [cmps finalizer exp config]
  {:c-files (get-final-and-tovalidate cmps finalizer config)
   :c1 {:name (:target finalizer)}
   :c2 {:name "validate"}
   :exp exp :dir (:dir config)})
