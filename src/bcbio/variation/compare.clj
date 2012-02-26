(ns bcbio.variation.compare
  "Generate comparisons between two sets of variant calls.
   Utilizes GATK walkers to generate detailed and summary statistics
   about two sets of calls:

   - Identify non-callable regions with CallableLociWalker
   - Combine variants from two samples
   - Use VariantEval to calculate overall concordance statistics
   - Provide output for concordant and discordant regions for
     detailed investigation"
  (:use [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template]]
        [bcbio.variation.stats :only [vcf-stats write-summary-table]]
        [bcbio.variation.report :only [concordance-report-metrics
                                       write-concordance-metrics
                                       write-scoring-table
                                       top-level-metrics]]
        [bcbio.variation.combine :only [combine-variants create-merged
                                        select-by-sample gatk-normalize]]
        [bcbio.variation.annotation :only [add-variant-annotations]]
        [bcbio.variation.filter :only [variant-filter pipeline-recalibration]]
        [bcbio.variation.phasing :only [is-haploid? compare-two-vcf-phased]]
        [clojure.math.combinatorics :only [combinations]]
        [clojure.java.io]
        [clojure.string :only [join]]
        [ordered.map :only [ordered-map]])
  (:require [fs.core :as fs]
            [clj-yaml.core :as yaml]
            [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]))

;; ## Variance assessment

(defn variant-comparison
  "Compare two variant files with GenotypeConcordance in VariantEval"
  [sample vcf1 vcf2 ref & {:keys [out-base interval-file]
                           :or {out-base nil iinterval-file nil}}]
  (let [file-info {:out-eval (str (itx/file-root (if (nil? out-base) vcf1 out-base)) ".eval")}
        args (concat
              ["-R" ref
               "--out" :out-eval
               "--eval" vcf1
               "--comp" vcf2
               "--sample" sample
               "--doNotUseAllStandardModules"
               "--evalModule" "CompOverlap"
               "--evalModule" "CountVariants"
               "--evalModule" "GenotypeConcordance"
               "--evalModule" "TiTvVariantEvaluator"
               "--evalModule" "ValidationReport"
               "--stratificationModule" "Sample"
               "--stratificationModule" "Filter"]
              (if-not (nil? interval-file) ["-L:bed" interval-file] []))]
    (broad/run-gatk "VariantEval" args file-info {:out [:out-eval]})
    (:out-eval file-info)))

(defn select-by-concordance
  "Variant comparison producing 3 files: concordant and both directions discordant"
  [sample call1 call2 ref & {:keys [out-dir interval-file]
                             :or {out-dir nil interval-file nil}}]
  (let [base-dir (if (nil? out-dir) (fs/parent (:file call1)) out-dir)]
    (if-not (fs/exists? base-dir)
      (fs/mkdirs base-dir))
    (doall
     (for [[c1 c2 cmp-type] [[call1 call2 "concordance"]
                             [call1 call2 "discordance"]
                             [call2 call1 "discordance"]]]
       (let [file-info {:out-vcf (str (fs/file base-dir
                                               (format "%s-%s%s-%s%s-%s.vcf"
                                                       sample (:name c1) (get c1 :mod "")
                                                       (:name c2) (get c2 :mod "") cmp-type)))}
             args (concat
                   ["-R" ref
                    "--sample_name" sample
                    "--variant" (:file c1)
                    (str "--" cmp-type) (:file c2)
                    "--out" :out-vcf]
                   (if-not (nil? interval-file) ["-L:bed" interval-file] []))]
         (broad/run-gatk "SelectVariants" args file-info {:out [:out-vcf]})
         (:out-vcf file-info))))))

;; ## Custom parsing and combinations
;; Utilizes GATK VariantContexts

(defn- vc-by-match-category
  "Lazy stream of VariantContexts categorized by concordant/discordant matching."
  [in-file]
  (letfn [(genotype-alleles [g]
            (vec (map #(.toString %) (:alleles g))))
          (is-concordant? [vc]
            (= (-> (map genotype-alleles (:genotypes vc))
                   set
                   count)
               1))]
    (for [vc (parse-vcf in-file)]
      [(if (is-concordant? vc) :concordant :discordant)
       (:vc vc)])))

(defn split-variants-by-match
  "Provide concordant and discordant variants for two variant files."
  [vcf1 vcf2 ref]
  (let [combo-file (combine-variants [vcf1 vcf2] ref)
        out-map {:concordant (itx/add-file-part combo-file "concordant")
                 :discordant (itx/add-file-part combo-file "discordant")}]
    (if-not (fs/exists? (:concordant out-map))
      (write-vcf-w-template combo-file out-map (vc-by-match-category combo-file)
                            ref))
    out-map))

;; ## Top-level
;; Process a directory of variant calls from multiple
;; sources, generating a summary of concordance plus detailed metrics
;; differences for tweaking filters.

(defn- get-summary-writer [config config-file ext]
  (if-not (nil? (get-in config [:dir :out]))
    (do
      (if-not (fs/exists? (get-in config [:dir :out]))
        (fs/mkdirs (get-in config :dir :out)))
      (writer (str (fs/file (get-in config [:dir :out])
                            (format "%s-%s"
                                    (itx/file-root (fs/base-name config-file)) ext)))))
    (writer System/out)))

(defn- prepare-vcf-calls
  "Prepare merged and annotated VCF files for an experiment."
  [exp config]
  (let [align-bams (map #(get % :align (:align exp)) (:calls exp))
        out-dir (get-in config [:dir :prep] (get-in config [:dir :out]))
        start-vcfs (map #(gatk-normalize % exp out-dir) (:calls exp))
        sample-vcfs (map #(select-by-sample (:sample exp) % (:ref exp)
                                            :intervals (:intervals exp)
                                            :out-dir out-dir)
                         start-vcfs)
        all-intervals (remove nil? (map :intervals (cons exp (:calls exp))))
        merged-vcfs (create-merged sample-vcfs align-bams (map #(get % :refcalls true) (:calls exp))
                                   (:ref exp) :out-dir (get-in config [:dir :out])
                                   :intervals all-intervals)
        ann-vcfs (map (fn [[v b c]] (if (get c :annotate false)
                                        (add-variant-annotations v b (:ref exp))
                                        v))
                      (map vector merged-vcfs align-bams (:calls exp)))
        filter-vcfs (map (fn [[v c]] (if-not (nil? (:filters c))
                                       (variant-filter v (:filters c) (:ref exp))
                                       v))
                         (map vector ann-vcfs (:calls exp)))]
    (map (fn [[c v]] (assoc c :file v))
         (map vector (:calls exp) filter-vcfs))))

(defn- compare-two-vcf-standard
  "Compare two standard VCF files based on the supplied configuration."
  [c1 c2 exp config]
  (let [c-files (select-by-concordance (:sample exp) c1 c2 (:ref exp)
                                       :out-dir (get-in config [:out :dir])
                                       :interval-file (:intervals exp))
        eval-file (variant-comparison (:sample exp) (:file c1) (:file c2)
                                      (:ref exp) :out-base (first c-files)
                                      :interval-file (:intervals exp))
        metrics (first (concordance-report-metrics (:sample exp) eval-file))]
    {:c-files c-files :metrics metrics :c1 c1 :c2 c2 :sample (:sample exp)}))

(defn- compare-two-vcf
  "Compare two VCF files, handling standard and haploid specific comparisons."
  [c1 c2 exp config]
  (let [phased-vcfs (group-by #(-> % :file is-haploid?) [c1 c2])]
    (if (get phased-vcfs true)
      (compare-two-vcf-phased (first (get phased-vcfs false))
                              (first (get phased-vcfs true))
                              exp config)
      (compare-two-vcf-standard c1 c2 exp config))))

(defn finalize-comparisons
  "Finalize all comparisons with finished initial pass data."
  [cmps exp config]
  (let [finalize-fns {"recal-filter" pipeline-recalibration}
        cmps-by-name (reduce (fn [m x] (assoc m [(-> x :c1 :name)
                                                 (-> x :c2 :name)] x))
                             (ordered-map)
                             cmps)]
    (letfn [(add-summary [x]
              (assoc x :summary (top-level-metrics x)))
            (update-w-finalizer [cur-cmps finalizer]
              "Update the current comparisons with a defined finalizer."
              (let [updated-cmp ((get finalize-fns (:method finalizer))
                                 (get cmps-by-name (:target finalizer))
                                 (get cmps-by-name (get finalizer :support (:target finalizer)))
                                 (:ref exp))]
                (assoc cur-cmps (:target finalizer)
                       (compare-two-vcf (:c1 updated-cmp) (:c2 updated-cmp) exp config))))]
      (map add-summary (vals (reduce update-w-finalizer cmps-by-name (:finalize exp)))))))

(defn variant-comparison-from-config
  "Perform comparison between variant calls using inputs from YAML config."
  [config-file]
  (let [config (-> config-file slurp yaml/parse-string)
        comparisons (flatten
                     (for [exp (:experiments config)]
                       (let [cmps (for [[c1 c2] (combinations (prepare-vcf-calls exp config) 2)]
                                    (compare-two-vcf c1 c2 exp config))]
                         (finalize-comparisons cmps exp config))))]
    (with-open [w (get-summary-writer config config-file "summary.txt")]
      (doseq [x comparisons]
        (.write w (format "* %s : %s vs %s\n" (:sample x)
                          (-> x :c1 :name) (-> x :c2 :name)))
        (write-scoring-table (:metrics x) w)
        (write-concordance-metrics (:summary x) w)
        (doseq [f (:c-files x)]
          (.write w (format "** %s\n" (fs/base-name f)))
          (write-summary-table (vcf-stats f) :wrtr w))))
    (with-open [w (get-summary-writer config config-file "summary.csv")]
      (doseq [[i x] (map-indexed vector (map :summary comparisons))]
        (if (= i 0)
          (.write w (format "%s\n" (join "," (map name (keys x))))))
        (.write w (format "%s\n" (join "," (for [v ( vals x)]
                                             (if (map? v) (:total v) v)))))))
    comparisons))

(defn -main [config-file]
  (variant-comparison-from-config config-file)
  nil)
