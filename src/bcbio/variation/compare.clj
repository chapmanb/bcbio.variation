(ns bcbio.variation.compare
  "Generate comparisons between two sets of variant calls.
   Utilizes GATK walkers to generate detailed and summary statistics
   about two sets of calls:

   - Identify non-callable regions with CallableLociWalker
   - Combine variants from two samples
   - Use VariantEval to calculate overall concordance statistics
   - Provide output for concordant and discordant regions for
     detailed investigation"
  (:use [clojure.java.io]
        [clojure.math.combinatorics :only [combinations]]
        [ordered.map :only [ordered-map]]
        [bcbio.align.reorder :only [reorder-bam]]
        [bcbio.variation.annotation :only [add-variant-annotations]]
        [bcbio.variation.callable :only [get-callable-bed]]
        [bcbio.variation.combine :only [combine-variants gatk-normalize]]
        [bcbio.variation.config :only [load-config do-transition]]
        [bcbio.variation.evaluate :only [calc-variant-eval-metrics]]
        [bcbio.variation.filter :only [variant-filter variant-format-filter
                                       pipeline-recalibration]]
        [bcbio.variation.filter.intervals :only [combine-multiple-intervals]]
        [bcbio.variation.multiple :only [prep-cmp-name-lookup pipeline-compare-multiple]]
        [bcbio.variation.multisample :only [compare-two-vcf-flexible
                                            multiple-samples?]]
        [bcbio.variation.recall :only [create-merged]]
        [bcbio.variation.structural :only [compare-sv-pipeline]]
        [bcbio.variation.validate :only [pipeline-validate]]
        [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template
                                               get-vcf-iterator]])
  (:require [clojure.string :as string]
            [clj-yaml.core :as yaml]
            [fs.core :as fs]
            [lonocloud.synthread :as ->]
            [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]
            [bcbio.variation.grade :as grade]
            [bcbio.variation.phasing :as phasing]
            [bcbio.variation.report :as report]))

;; ## Variance assessment

(defn select-by-concordance
  "Variant comparison producing 3 files: concordant and both directions discordant"
  [sample call1 call2 ref & {:keys [out-dir intervals]}]
  (let [base-dir (if (nil? out-dir) (fs/parent (:file call1)) out-dir)
        ready-intervals (remove nil? (flatten [intervals (:intervals call1)
                                               (:intervals call2)]))]
    (if-not (fs/exists? base-dir)
      (fs/mkdirs base-dir))
    (doall
     (for [[c1 c2 cmp-type] [[call1 call2 "concordance"]
                             [call1 call2 "discordance"]
                             [call2 call1 "discordance"]]]
       (let [file-info {:out-vcf (str (fs/file base-dir
                                               (format "%s-%s-%s-%s.vcf"
                                                       sample (:name c1) (:name c2) cmp-type)))}
             args (concat
                   ["-R" ref
                    "--sample_name" sample
                    "--variant" (:file c1)
                    (str "--" cmp-type) (:file c2)
                    "--out" :out-vcf]
                   (broad/gatk-cl-intersect-intervals ready-intervals ref))]
         (broad/run-gatk "SelectVariants" args file-info {:out [:out-vcf]})
         (:out-vcf file-info))))))

;; ## Custom parsing and combinations
;; Utilizes GATK VariantContexts

(defn- vc-by-match-category
  "Lazy stream of VariantContexts categorized by concordant/discordant matching."
  [vcf-iter]
  (letfn [(genotype-alleles [g]
            (vec (map #(.toString %) (:alleles g))))
          (is-concordant? [vc]
            (= (-> (map genotype-alleles (:genotypes vc))
                   set
                   count)
               1))]
    (for [vc (parse-vcf vcf-iter)]
      [(if (is-concordant? vc) :concordant :discordant)
       (:vc vc)])))

(defn split-variants-by-match
  "Provide concordant and discordant variants for two variant files."
  [vcf1 vcf2 ref]
  (let [combo-file (combine-variants [vcf1 vcf2] ref)
        out-map {:concordant (itx/add-file-part combo-file "concordant")
                 :discordant (itx/add-file-part combo-file "discordant")}]
    (if-not (fs/exists? (:concordant out-map))
      (with-open [combo-vcf-iter (get-vcf-iterator combo-file ref)]
        (write-vcf-w-template combo-file out-map (vc-by-match-category combo-vcf-iter)
                              ref)))
    out-map))

;; ## Pipeline
;; Process a directory of variant calls from multiple
;; sources, generating a summary of concordance plus detailed metrics
;; differences for tweaking filters.

(defn- prepare-input-bams
  "Retrieve BAM files associated with alignments, normalizing if needed."
  [exp out-dir]
  (let [call-bams (map (fn [c] [(get c :align (:align exp)) c]) (:calls exp))]
    (map (fn [[b c]] (when-not (nil? b)
                     (reorder-bam b (:ref exp) c exp :out-dir out-dir)))
         call-bams)))

(defn- prepare-vcf-calls
  "Prepare merged and annotated VCF files for an experiment."
  [exp config]
  (let [out-dir (get-in config [:dir :prep] (get-in config [:dir :out]))
        transition (partial do-transition config)
        align-bams (prepare-input-bams exp out-dir)
        all-intervals (remove nil? (map :intervals (cons exp (:calls exp))))
        start-vcfs (vec (map #(gatk-normalize % exp all-intervals out-dir transition)
                             (:calls exp)))
        _ (transition :combine "Creating merged VCF files for all comparisons")
        merged-vcfs (create-merged (map :file start-vcfs) align-bams exp
                                   :out-dir out-dir
                                   :intervals all-intervals)
        _ (transition :annotate "Annotate VCFs with metrics")
        ann-vcfs (map (fn [[v b c]]
                        (add-variant-annotations v b (:ref exp) c :out-dir out-dir
                                                 :intervals all-intervals))
                      (map vector merged-vcfs align-bams (:calls exp)))
        _ (transition :filter "Post annotation filtering")
        filter-vcfs (map (fn [[v c]]
                           (cond
                            (:filters c) (variant-filter v (:filters c) (:ref exp))
                            (:format-filters c) (variant-format-filter v (:format-filters c)
                                                                       (:ref exp))
                            :else v))
                         (map vector ann-vcfs (:calls exp)))]
    (map (fn [[c v b]] (-> c
                           (assoc :file v)
                           (assoc :align b)))
         (map vector (:calls exp) filter-vcfs align-bams))))

(defn- compare-two-vcf-standard
  "Compare two standard VCF files based on the supplied configuration."
  [c1 c2 exp config]
  (letfn [(callable-intervals [exp c1 c2]
            (let [out-dir (get-in config [:dir :prep] (get-in config [:dir :out]))
                  align-bams (remove nil? (map :align [c1 c2]))]
              (when (and (:intervals exp) (seq align-bams))
                (combine-multiple-intervals (:intervals exp) align-bams (:ref exp)
                                            :out-dir out-dir :name (:sample exp)))))
          (discordant-name [x]
            (format "%s-discordant" (:name x)))
          (zipmap-ordered [xs1 xs2]
            (apply ordered-map (interleave xs1 xs2)))]
    (let [c-files (select-by-concordance (:sample exp) c1 c2 (:ref exp)
                                         :out-dir (get-in config [:dir :out])
                                         :intervals (:intervals exp))
          eval (calc-variant-eval-metrics (:sample exp) (:file c1) (:file c2) (:ref exp)
                                          :out-base (first c-files)
                                          :intervals (:intervals exp))
          c-eval (calc-variant-eval-metrics (:sample exp) (:file c1) (:file c2) (:ref exp)
                                            :out-base (itx/add-file-part (first c-files) "callable")
                                            :intervals (callable-intervals exp c1 c2))]
      {:c-files (zipmap-ordered (map keyword
                                     ["concordant" (discordant-name c1) (discordant-name c2)])
                                c-files)
       :c1 c1 :c2 c2 :exp exp :dir (config :dir)
       :metrics (report/concordance-report-metrics (:sample exp) eval)
       :callable-metrics (report/concordance-report-metrics (:sample exp) c-eval)})))

(defn compare-two-vcf
  "Compare two VCF files, handling standard and haploid specific comparisons."
  [c1 c2 exp config]
  (do-transition config :compare (format "Comparing VCFs: %s vs %s" (:name c1) (:name c2)))
  (let [[c1 c2 sv-cmp] (if-not (:mod c1)
                         (compare-sv-pipeline c1 c2 exp config)
                         [c1 c2 {}])
        phased-vcfs (group-by #(-> % :file (phasing/is-haploid? (:ref exp))) [c1 c2])
        out-cmp (cond
                 (get phased-vcfs true) (phasing/compare-two-vcf-phased phased-vcfs exp config)
                 (multiple-samples? (:file c1)) (compare-two-vcf-flexible c1 c2 exp config)
                 :else (compare-two-vcf-standard c1 c2 exp config))
        grade-cmp (if (grade/is-grade-cmp? exp)
                    (grade/annotate-discordant out-cmp)
                    out-cmp)]
    (assoc grade-cmp :c-files (reduce (fn [coll [k v]] (assoc coll k v))
                                      (:c-files grade-cmp) sv-cmp))))

;; ## Customizable finalizer comparisons

(defmulti run-finalizer
  "Run a post-pairwise comparison function, returning updated comparison details,"
  (fn [cmps finalizer exp config] (-> finalizer :method keyword)))

(defmethod run-finalizer :recal-filter
  [& args]
  (apply pipeline-recalibration args))

(defmethod run-finalizer :multiple
  [& args]
  (apply pipeline-compare-multiple args))

(defmethod run-finalizer :validate
  [& args]
  (apply pipeline-validate args))

(defn finalize-comparisons
  "Finalize all comparisons with finished initial pass data."
  [cmps exp config]
  (letfn [(add-summary [x]
            (-> x
                (assoc :exp exp)
                (->/as cur-cmp
                  (assoc :summary (report/top-level-metrics cur-cmp)))
                (->/when (grade/is-grade-cmp? exp)
                  (->/as cmp-w-summary
                    (assoc :grade-breakdown (grade/prep-grade-breakdown cmp-w-summary))))))
          (update-w-finalizer [cur-cmps finalizer]
            "Update the current comparisons with a defined finalizer."
            (do-transition config :finalize
                           (format "Finalize %s: %s" (:method finalizer)
                                   (let [t (:target finalizer)]
                                     (if (coll? t) (string/join ", " t) t))))
            (let [updated-cmp (run-finalizer cur-cmps finalizer exp config)]
              (assoc cur-cmps (map #(get-in updated-cmp [% :name]) [:c1 :c2])
                     (if-not (:re-compare updated-cmp) updated-cmp
                             (compare-two-vcf (:c1 updated-cmp) (:c2 updated-cmp) exp config)))))]
    (->> (reduce update-w-finalizer
                 (prep-cmp-name-lookup cmps) (:finalize exp))
         vals
         (map add-summary))))

;; ## Top-level

(defn- get-summary-writer [config config-file ext]
  (if-not (nil? (get-in config [:dir :out]))
    (do
      (if-not (fs/exists? (get-in config [:dir :out]))
        (fs/mkdirs (get-in config :dir :out)))
      (writer (str (fs/file (get-in config [:dir :out])
                            (format "%s-%s"
                                    (itx/file-root (fs/base-name config-file)) ext)))))
    (writer System/out)))

(defn variant-comparison-from-config
  "Perform comparison between variant calls using inputs from YAML config."
  [config-file]
  (let [config (load-config config-file)
        comparisons (flatten
                     (for [exp (:experiments config)]
                       (let [cmps (for [[c1 c2] (combinations (prepare-vcf-calls exp config) 2)]
                                    (compare-two-vcf c1 c2 exp config))]
                         (finalize-comparisons cmps exp config))))
        grading-file (str (fs/file (get-in config [:dir :out])
                                   (format "%s-grading.yaml" (itx/file-root (fs/base-name config-file)))))]
    (do-transition config :summary "Summarize comparisons")
    (with-open [w (get-summary-writer config config-file "summary.txt")]
      (report/write-summary-txt w comparisons))
    (with-open [w (get-summary-writer config config-file "files.csv")]
      (report/write-files-csv w comparisons config))
    (with-open [w (get-summary-writer config config-file "summary.csv")]
      (report/write-summary-csv w comparisons))
    (when-let [bdowns (seq (remove nil? (map :grade-breakdown comparisons)))]
      (spit grading-file (yaml/generate-string bdowns)))
    (do-transition config :finished "Finished")
    comparisons))

(defn -main [config-file]
  (try
    (variant-comparison-from-config config-file)
    (catch Throwable t
      (.printStackTrace t)
      (System/exit -1))
    (finally
      (shutdown-agents)
      (System/exit 0))))
