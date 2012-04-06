(ns bcbio.variation.compare
  "Generate comparisons between two sets of variant calls.
   Utilizes GATK walkers to generate detailed and summary statistics
   about two sets of calls:

   - Identify non-callable regions with CallableLociWalker
   - Combine variants from two samples
   - Use VariantEval to calculate overall concordance statistics
   - Provide output for concordant and discordant regions for
     detailed investigation"
  (:use [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template
                                               get-vcf-source]]
        [bcbio.variation.metrics :only [vcf-stats write-summary-table]]
        [bcbio.variation.report :only [concordance-report-metrics
                                       write-concordance-metrics
                                       write-scoring-table
                                       top-level-metrics
                                       write-classification-metrics]]
        [bcbio.variation.combine :only [combine-variants create-merged
                                        gatk-normalize gatk-cl-intersect-intervals]]
        [bcbio.variation.annotation :only [add-variant-annotations]]
        [bcbio.variation.filter :only [variant-filter pipeline-recalibration]]
        [bcbio.variation.phasing :only [is-haploid? compare-two-vcf-phased]]
        [bcbio.variation.callable :only [get-callable-bed]]
        [bcbio.variation.multiple :only [prep-cmp-name-lookup pipeline-compare-multiple]]
        [bcbio.align.reorder :only [reorder-bam]]
        [ordered.map :only [ordered-map]]
        [clojure.math.combinatorics :only [combinations]]
        [clojure.java.io])
  (:require [clojure.string :as string]
            [clojure.data.csv :as csv]
            [fs.core :as fs]
            [clj-yaml.core :as yaml]
            [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]))

;; ## Variance assessment

(defn calc-variant-eval-metrics
  "Compare two variant files with GenotypeConcordance in VariantEval"
  [sample vcf1 vcf2 ref & {:keys [out-base intervals]}]
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
              (gatk-cl-intersect-intervals intervals))]
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
                                               (format "%s-%s-%s-%s.vcf"
                                                       sample (:name c1) (:name c2) cmp-type)))}
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
  [vcf-source]
  (letfn [(genotype-alleles [g]
            (vec (map #(.toString %) (:alleles g))))
          (is-concordant? [vc]
            (= (-> (map genotype-alleles (:genotypes vc))
                   set
                   count)
               1))]
    (for [vc (parse-vcf vcf-source)]
      [(if (is-concordant? vc) :concordant :discordant)
       (:vc vc)])))

(defn split-variants-by-match
  "Provide concordant and discordant variants for two variant files."
  [vcf1 vcf2 ref]
  (let [combo-file (combine-variants [vcf1 vcf2] ref)
        out-map {:concordant (itx/add-file-part combo-file "concordant")
                 :discordant (itx/add-file-part combo-file "discordant")}]
    (if-not (fs/exists? (:concordant out-map))
      (with-open [combo-vcf-s (get-vcf-source combo-file ref)]
        (write-vcf-w-template combo-file out-map (vc-by-match-category combo-vcf-s)
                              ref)))
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
        align-bams (prepare-input-bams exp out-dir)
        start-vcfs (map #(gatk-normalize % exp out-dir) (:calls exp))
        all-intervals (remove nil? (map :intervals (cons exp (:calls exp))))
        merged-vcfs (create-merged (map :file start-vcfs)
                                   align-bams
                                   (map #(get % :refcalls true) (:calls exp))
                                   (:ref exp) :out-dir out-dir 
                                   :intervals all-intervals)
        ann-vcfs (map (fn [[v b c]] (if (get c :annotate false)
                                        (add-variant-annotations v b (:ref exp) :out-dir out-dir)
                                        v))
                      (map vector merged-vcfs align-bams (:calls exp)))
        filter-vcfs (map (fn [[v c]] (if-not (nil? (:filters c))
                                       (variant-filter v (:filters c) (:ref exp))
                                       v))
                         (map vector ann-vcfs (:calls exp)))]
    (map (fn [[c v b]] (-> c
                           (assoc :file v)
                           (assoc :align b)))
         (map vector (:calls exp) filter-vcfs align-bams))))

(defn- compare-two-vcf-standard
  "Compare two standard VCF files based on the supplied configuration."
  [c1 c2 exp config]
  (letfn [(callable-intervals [exp c1 c2]
            (let [out-dir (get-in config [:dir :prep] (get-in config [:dir :out]))]
              (remove nil? (cons (:intervals exp)
                                 (map #(when-not (nil? (:align %))
                                         (get-callable-bed (:align %) (:ref exp)
                                                           :out-dir out-dir))
                                      [c1 c2])))))
          (discordant-name [x]
            (format "%s-discordant" (:name x)))
          (zipmap-ordered [xs1 xs2]
            (apply ordered-map (interleave xs1 xs2)))]
    (let [c-files (select-by-concordance (:sample exp) c1 c2 (:ref exp)
                                         :out-dir (get-in config [:dir :out])
                                         :interval-file (:intervals exp))
          eval (calc-variant-eval-metrics (:sample exp) (:file c1) (:file c2) (:ref exp)
                                          :out-base (first c-files)
                                          :intervals (:intervals exp))
          c-eval (calc-variant-eval-metrics (:sample exp) (:file c1) (:file c2) (:ref exp)
                                            :out-base (itx/add-file-part (first c-files) "callable")
                                            :intervals (callable-intervals exp c1 c2))]
      {:c-files (zipmap-ordered ["concordant" (discordant-name c1) (discordant-name c2)]
                                c-files)
       :c1 c1 :c2 c2 :exp exp :dir (config :dir)
       :metrics (first (concordance-report-metrics (:sample exp) eval))
       :callable-metrics (first (concordance-report-metrics (:sample exp) c-eval))})))

(defn- compare-two-vcf
  "Compare two VCF files, handling standard and haploid specific comparisons."
  [c1 c2 exp config]
  (let [phased-vcfs (group-by #(-> % :file (is-haploid? (:ref exp))) [c1 c2])]
    (if (get phased-vcfs true)
      (compare-two-vcf-phased (first (get phased-vcfs false))
                              (first (get phased-vcfs true))
                              exp config)
      (compare-two-vcf-standard c1 c2 exp config))))

(defn finalize-comparisons
  "Finalize all comparisons with finished initial pass data."
  [cmps exp config]
  (let [finalize-fns {"recal-filter" pipeline-recalibration
                      "multiple" pipeline-compare-multiple}
        cmps-by-name (prep-cmp-name-lookup cmps)]
    (letfn [(add-summary [x]
              (-> x
                  (assoc :exp exp)
                  (#(assoc % :summary (top-level-metrics %)))))
            (update-w-finalizer [cur-cmps finalizer]
              "Update the current comparisons with a defined finalizer."
              (let [updated-cmp ((get finalize-fns (:method finalizer))
                                 cmps-by-name finalizer exp config)]
                (assoc cur-cmps (map #(get-in updated-cmp [% :name]) [:c1 :c2])
                       (if-not (:re-compare updated-cmp) updated-cmp
                               (compare-two-vcf (:c1 updated-cmp) (:c2 updated-cmp) exp config)))))]
      (map add-summary (vals (reduce update-w-finalizer cmps-by-name (:finalize exp)))))))

(defn load-config
  "Load configuration file, handling conversion of relative to absolute paths."
  [config-file]
  (let [config (-> config-file slurp yaml/parse-string)
        base-dir (fs/file (get-in config [:dir :base] "."))
        to-process #{[:dir :out] [:dir :prep]
                     [:experiments :ref] [:experiments :intervals]
                     [:experiments :align] [:experiments :calls :file]
                     [:experiments :calls :align]}]
    (letfn [(make-absolute [x]
              (if (.isAbsolute (file x))
                x
                (str (fs/file base-dir x))))
            (maybe-process [val path]
              (if (contains? to-process path)
                (if (seq? val)
                  (map make-absolute val)
                  (make-absolute val))
                val))
            (update-tree [config path]
              (cond (map? config)
                    (reduce (fn [item [k v]]
                              (assoc item k (cond
                                             (map? v) (update-tree v (conj path k))
                                             (seq? v) (map #(update-tree % (conj path k)) v)
                                             :else (maybe-process v (conj path k)))))
                            config
                            (vec config))
                    (contains? to-process path) (maybe-process config path)
                    :else config))]
      (update-tree config []))))

(defn variant-comparison-from-config
  "Perform comparison between variant calls using inputs from YAML config."
  [config-file]
  (let [config (load-config config-file)
        comparisons (flatten
                     (for [exp (:experiments config)]
                       (let [cmps (for [[c1 c2] (combinations (prepare-vcf-calls exp config) 2)]
                                    (compare-two-vcf c1 c2 exp config))]
                         (finalize-comparisons cmps exp config))))]
    (with-open [w (get-summary-writer config config-file "summary.txt")
                w2 (get-summary-writer config config-file "files.csv")]
      (csv/write-csv w2 [["call1" "call2" "type" "fname"]])
      (doseq [x comparisons]
        (.write w (format "* %s : %s vs %s\n" (-> x :exp :sample)
                          (-> x :c1 :name) (-> x :c2 :name)))
        (write-scoring-table (:metrics x) w)
        (write-concordance-metrics (:summary x) w)
        (when (get-in x [:c1 :mod])
          (write-classification-metrics x w))
        (doseq [[k f] (:c-files x)]
          (.write w (format "** %s\n" (name k)))
          (csv/write-csv w2 [[(get-in x [:c1 :name]) (get-in x [:c2 :name]) (name k)
                              (string/replace-first f (str (get-in config [:dir :out]) "/") "")]])
          (write-summary-table (vcf-stats f (get-in x [:exp :ref])) :wrtr w))))
    (with-open [w (get-summary-writer config config-file "summary.csv")]
      (doseq [[i x] (map-indexed vector (map :summary comparisons))]
        (when (= i 0)
          (.write w (format "%s\n" (string/join "," (map name (keys x))))))
        (.write w (format "%s\n" (string/join "," (for [v (vals x)]
                                                    (if (map? v) (:total v) v)))))))
    comparisons))

(defn -main [config-file]
  (variant-comparison-from-config config-file)
  nil)
