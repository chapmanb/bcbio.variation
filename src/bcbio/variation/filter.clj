(ns bcbio.variation.filter
  "Filter variant calls according to supplied criteria."
  (:import [org.broadinstitute.sting.utils.variantcontext
            VariantContextBuilder])
  (:use [clojure.string :only [split]]
        [bcbio.variation.multiple :only [multiple-overlap-analysis]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template
                                               get-vcf-source]])
  (:require [clojure.string :as string]
            [incanter.stats :as stats]
            [bcbio.run.broad :as broad]
            [bcbio.run.itx :as itx]))

(defn jexl-from-config [jexl-filters]
  "Retrieve GATK JEXL commandline expressions from filters."
  (letfn [(jexl-args [x]
            ["--filterName" (str (first (split x #"\s+")) "Filter")
             "--filterExpression" x])]
    (flatten (map jexl-args jexl-filters))))

(defn variant-filter
  "Perform hard variant filtering with supplied JEXL expression criteria."
  [in-vcf jexl-filters ref]
  (let [file-info {:out-vcf (itx/add-file-part in-vcf "filter")}
        args (concat ["-R" ref
                      "--variant" in-vcf
                      "-o" :out-vcf
                      "-l" "ERROR"]
                      (jexl-from-config jexl-filters))]
    (broad/run-gatk "VariantFiltration" args file-info {:out [:out-vcf]})
    (:out-vcf file-info)))

(defn- variant-recalibration
  "Perform the variant recalibration step with input training VCF files.
  training-vcfs is a list of `{:file vcf-file :name name-to-use :prior probability}`"
  [in-vcf training-vcfs annotations ref & {:keys [lenient]}]
  (let [base-out (itx/file-root in-vcf)
        file-info {:out-recal (str base-out ".recal")
                   :out-tranch (str base-out ".tranches")
                   :out-r (str base-out "-recalplots.R")}
        args (concat ["-R" ref
                      "-input" in-vcf
                      "-recalFile" :out-recal
                      "-tranchesFile" :out-tranch
                      "-rscriptFile" :out-r
                      "--mode" "BOTH"]
                     (if lenient
                       ["--percentBadVariants" "0.05"
                        "--maxGaussians" "4"]
                       ["--percentBadVariants" "0.03"
                        "--maxGaussians" "10"])
                     (flatten (map (fn [x] ["-an" x]) annotations))
                     (flatten (map (fn [x] [(str "-resource:" (:name x)
                                                 ",known=true"
                                                 ",training=true"
                                                 ",truth=" (:truth x)
                                                 ",prior=" (:prior x)
                                                 ",bad=" (:bad x))
                                            (:file x)])
                                   training-vcfs)))]
    (broad/run-gatk "VariantRecalibrator" args file-info {:out [:out-recal :out-tranch]})
    file-info))

(defn- apply-recalibration
  "Apply variant recalibration to input VCF."
  [in-vcf recal-files ref]
  (let [file-info {:out-vcf (itx/add-file-part in-vcf "recalfilter")}
        args ["-R" ref
              "-input" in-vcf
              "--ts_filter_level" "99.0"
              "--mode" "BOTH"
              "-tranchesFile" (:out-tranch recal-files)
              "-recalFile" (:out-recal recal-files)
              "-o" :out-vcf]]
    (broad/run-gatk "ApplyRecalibration" args file-info {:out [:out-vcf]})
    (:out-vcf file-info)))

(defn variant-recal-filter
  "Perform filtration using variant recalibration based on known variations.
  Training-vcfs is a list of true training sites along with associated
  probability and name."
  [in-vcf training-vcfs annotations ref & {:keys [lenient]}]
  (let [recal-files (variant-recalibration in-vcf training-vcfs annotations ref :lenient lenient)]
    (apply-recalibration in-vcf recal-files ref)))

(defn remove-cur-filters
  "Remove any filter information in the supplied file."
  [in-vcf ref]
  (letfn [(remove-vc-filter [vc]
            [:out (-> (VariantContextBuilder. (:vc vc))
                      (.passFilters)
                      (.make))])]
    (let [out-file (itx/add-file-part in-vcf "nofilter")]
      (with-open [vcf-source (get-vcf-source in-vcf ref)]
        (write-vcf-w-template in-vcf {:out out-file}
                              (map remove-vc-filter (parse-vcf vcf-source))
                              ref))
      out-file)))

(defn- get-train-info
  "Retrieve training information for GATK recalibration:
   - No support specified: use the target comparison
   - Support specified and a specific comparison pair
   - Support specified as a single target: use target versus all comparison"
  [cmps-by-name support config]
  (let [support-vcfs (if (coll? support)
                       (take 2 (-> cmps-by-name (get support) :c-files vals))
                       (let [x (multiple-overlap-analysis cmps-by-name config support)]
                         [(:true-positives x) (:false-positives x)]))]
      [{:file (first support-vcfs)
        :name "concordant"
        :truth "true"
        :bad "false"
        :prior 10.0}
       {:file (second support-vcfs)
        :name "discordant"
        :truth "false"
        :bad "true"
        :prior 10.0}]))

(defn pipeline-recalibration
  "Perform variant recalibration and filtration as part of processing pipeline."
  [cmps-by-name finalizer exp config]
  (let [init-target (get cmps-by-name (:target finalizer))
        all-params (let [x (:params finalizer)] (if (map? x) [x] x))]
    (reduce (fn [target [params fkey]]
              (let [in-vcf (remove-cur-filters (-> target fkey :file) (:ref exp))
                    train-info (get-train-info cmps-by-name
                                               (get params :support (:target finalizer))
                                               config)]
                (-> target
                    (assoc-in [fkey :file]
                              (-> in-vcf
                                  (#(if-let [anns (:annotations params)]
                                      (variant-recal-filter % train-info
                                                            anns (:ref exp)
                                                            :lenient (:lenient params))
                                      %))
                                  (#(if-let [hard-filters (:filters params)]
                                      (variant-filter % hard-filters (:ref exp))
                                      %))))
                    (#(assoc-in % [fkey :name] (format "%s-%s" (get-in % [fkey :name]) "recal")))
                    (assoc-in [fkey :mod] "recal")
                    (assoc :re-compare true))))
            init-target (map vector all-params [:c1 :c2]))))

;; ## Normalized attribute access

(defmulti get-vc-attr
  "Generalized retrieval of attributes from variant with a single genotype."
  (fn [vc attr] attr))

(defmethod get-vc-attr "AD"
  [vc attr]
  "AD: Allelic depth for ref and alt alleles. Converted to percent
   deviation from expected for haploid/diploid calls."
  (let [g (-> vc :genotypes first)
        ads (map #(Integer/parseInt %) (string/split (get-in g [:attributes attr]) #","))
        alleles (cons (:ref-allele vc) (:alt-alleles vc))
        ref-count (first ads)
        allele-count (apply + (map #(nth ads (.indexOf alleles %)) (set (:alleles g))))]
    (when-let [e-pct (get {"HOM_VAR" 1.0 "HET" 0.5 "HOM_REF" 0.0} (:type g))]
      (Math/abs (- e-pct (/ allele-count (+ allele-count ref-count)))))))

(defmethod get-vc-attr "QUAL"
  [vc attr]
  (:qual vc))

(defmethod get-vc-attr :default
  [vc attr]
  (let [x (get-in vc [:attributes attr])]
    (try (Float/parseFloat x)
         (catch java.lang.NumberFormatException _ x))))

(defn get-vc-attrs
  "Retrieve attributes from variants independent of location."
  [vc attrs]
  {:pre [(= 1 (count (:genotypes vc)))
         (contains? #{1 2} (-> vc :genotypes first :alleles count))]}
  (zipmap attrs (map (partial get-vc-attr vc) attrs)))

(defn get-vc-attr-ranges
  "Retrieve first/third quantile ranges of attributes for min/max normalization."
  [attrs in-vcf ref]
  (letfn [(get-quartiles [[k v]]
            [k (stats/quantile v :probs [0.25 0.75])])]
    (with-open [vcf-s (get-vcf-source in-vcf ref)]
      (->> (reduce (fn [coll vc]
                    (reduce (fn [icoll [k v]]
                              (assoc icoll k (cons v (get icoll k))))
                            coll (get-vc-attrs vc attrs)))
                  (zipmap attrs (repeat [])) (parse-vcf vcf-s))
           (map get-quartiles)
           (into {})))))

(defn get-vc-attrs-normalized
  "Min-max Normalized attributes for each variant context in an input file."
  [attrs in-vcf ref]
  (letfn [(min-max-norm [x [minv maxv]]
            (let [trunc-score-max (if (< x maxv) x maxv)
                  trunc-score (if (> trunc-score-max minv) trunc-score-max minv)]
              (/ (- trunc-score minv) (- maxv minv))))
          (min-max-norm-ranges [mm-ranges [k v]]
            [k (min-max-norm v (get mm-ranges k))])]
    (let [mm-ranges (get-vc-attr-ranges attrs in-vcf ref)]
      (fn [vc]
        (->> (get-vc-attrs vc attrs)
             (map (partial min-max-norm-ranges mm-ranges))
             (into {}))))))
