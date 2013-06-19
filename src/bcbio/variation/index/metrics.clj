(ns bcbio.variation.index.metrics
  "Pre-index a variant file for quick retrieval of associated metrics."
  (:use [ordered.map :only [ordered-map]]
        [bcbio.variation.metrics :only [passes-filter?]]
        [bcbio.variation.filter.attr :only [get-vc-attrs]]
        [bcbio.variation.index.subsample :only [subsample-by-cluster]]
        [bcbio.variation.variantcontext :only [get-vcf-header get-vcf-iterator parse-vcf]]
        [bcbio.variation.web.db :only [get-sqlite-db]])
  (:require [clojure.string :as string]
            [clojure.java.jdbc :as sql]
            [me.raynes.fs :as fs]
            [bcbio.run.itx :as itx]))

(def ^{:doc "Metrics to expose, ranked in order of priority with default min/max values."}
  expose-metrics
  (ordered-map "QUAL" {:range [0.0 2000.0]
                       :desc "Variant quality score, phred-scaled"}
               "DP" {:range [0.0 250.0]
                     :desc "Read depth after filtering of low quality reads"}
               "MQ" {:range [25.0 75.0]
                     :desc "Mapping quality"}
               "QD" {:range [0.0 50.0]
                     :desc "Variant confidence by depth"}
               "FS" {:range [0.0 100.0]
                     :y-scale {:type :log}
                     :desc "Phred-scaled p-value using Fisher's exact test to detect strand bias"}
               "ReadPosEndDist" {:range [0.0 50.0]
                                 :desc "Mean distance from either end of read"}
               "MFE" {:range [-10.0 0.0]
                      :y-scale {:type :log}
                      :desc (str "Minimum Free Energy of secondary structure near variant. "
                                 "Larger negative values are more problematic.")}
               "Entropy" {:range [1.0 4.0]
                          :desc (str "Shannon entropy of variant flanking regions. "
                                     "Low values indicate repetitive sequence.")}
               "AD" {:range [0.0 1.0]
                     :y-scale {:type :log}
                     :desc "Deviation from expected allele balance for ref/alt alleles"}
               "PL" {:range [-100.0 0]
                     :desc "Normalized, phred-scaled likelihoods for alternative genotype"}
               "HaplotypeScore" {:range [0.0 50.0]
                                 :desc "Consistency of the site with at most two segregating haplotypes"}
               ;; Validation
               "GradeCat" {:x-scale {:type :category}
                           :desc "Validation category, differentiating discordant types"}
               ;; Methylation metrics from BisSNP
               "CS" {:x-scale {:type :category}
                     :desc "Strand of cytosine relative to reference genome"}
               "Context" {:x-scale {:type :category}
                          :desc "Cytosine context: homozygous or heterozygous CG sites"}
               "CM" {:range [0.0 50.0]
                     :y-scale {:type :log}
                     :desc "Number of unconverted, methylated, cytosines"}
               "CU" {:range [0.0 1.0]
                     :y-scale {:type :log}
                     :desc "Percentage of methylated bases at a position."}))

(def ^{:doc "Default metrics that are always available." :private true}
  default-metrics
  [{:id "QUAL" :type :float}])

(defmulti available-metrics
  (fn [in-file] in-file))

(defmethod available-metrics nil
  ^{:doc "Retrieve all available default metrics without file information"}
  [_]
  (map (fn [[k v]] (assoc v :id k)) expose-metrics))

(defmethod available-metrics :default
  ^{:doc "Retrieve metrics available for variant input file."}
  [vcf-file]
  (letfn [(convert-header [line]
            {:id (.getID line)
             :type (case (.name (.getType line))
                     "Integer" :float
                     "Float" :float
                     "String" :text
                     "Character" :text
                     :else nil)})
          (add-base-info [x]
            (merge x (get expose-metrics (:id x))))]
    (let [metrics-order (reduce (fn [coll [i x]] (assoc coll x i))
                                {} (map-indexed vector (keys expose-metrics)))]
      (->> (get-vcf-header vcf-file)
           .getMetaDataInInputOrder
           (filter #(contains? #{"INFO" "FORMAT"} (.getKey %)))
           (filter #(contains? expose-metrics (.getID %)))
           (group-by #(.getID %))
           vals
           (map first)
           (map convert-header)
           (concat default-metrics)
           (map add-base-info)
           (sort-by #(get metrics-order (:id %)))))))

(def ^{:doc "Common columns for variant metrics table."
       :private true}
  shared-metrics-cols [[:contig :text]
                       [:start :integer]
                       [:refallele :text]
                       [:issubsample :integer]])

(defn- create-metrics-tables
  "Create table to represent variant metrics"
  [metrics]
  (apply sql/create-table (concat [:metrics] shared-metrics-cols
                                  (map (fn [x] [(:id x) (:type x)]) metrics))))

(defn- index-needs-update?
  "Check if an index file has up to date columns"
  [index-file metrics]
  (let [want-cols (set (concat (map first shared-metrics-cols)
                               (map #(keyword (string/lower-case (:id %))) metrics)))]
    (sql/with-connection (get-sqlite-db index-file)
      (sql/with-query-results rows
        ["SELECT * from metrics LIMIT 1"]
        (not= want-cols (-> rows first keys set))))))

(defn- num-variants
  [index-file]
  (sql/with-connection (get-sqlite-db index-file)
    (sql/with-query-results rows
      ["SELECT count(*) from metrics"]
      (-> rows first vals first))))

(defn- all-metrics-as-subsample
  [index-file]
  (sql/with-connection (get-sqlite-db index-file)
    (sql/update-values :metrics ["start > -1"] {:issubsample 1})))

(declare get-raw-metrics)
(defn get-raw-metrics-linear
  "Retrieve numerical raw metrics for filtering"
  [in-file ref-file]
  (let [metrics (->> (available-metrics in-file)
                     (filter #(= :linear (get-in % [:x-scale :type] :linear)))
                     (map :id))]
    (get-raw-metrics in-file ref-file :metrics metrics)))

(defn- subsample-metrics
  "Identify a subsample of records to use in visualization"
  [index-file in-file ref-file params]
  (let [sub-ids (subsample-by-cluster (get-raw-metrics-linear in-file ref-file) params)]
    (sql/with-connection (get-sqlite-db index-file)
      (sql/transaction
       (doseq [xid sub-ids]
         (sql/update-values :metrics
                            (cons "contig=? AND start=? and refallele=?" (vec xid))
                            {:issubsample 1}))))))

(defn index-variant-file
  "Pre-index a variant file with associated metrics."
  [in-file ref-file & {:keys [re-index? subsample-params]}]
  (let [batch-size 10000
        metrics (available-metrics in-file)
        index-file (str (itx/file-root in-file) "-metrics.db")]
    (when (or re-index?
              (itx/needs-run? index-file)
              (index-needs-update? index-file metrics))
      (itx/with-tx-file [tx-index index-file]
        (sql/with-connection (get-sqlite-db tx-index :create true)
          (sql/transaction
           (create-metrics-tables metrics))
          (with-open [vcf-iter (get-vcf-iterator in-file ref-file)]
            (doseq [vcs (partition-all batch-size (filter passes-filter? (parse-vcf vcf-iter)))]
              (sql/transaction
               (doseq [vc vcs]
                 (sql/insert-record :metrics
                                    (-> (get-vc-attrs vc (map :id metrics) {})
                                        (assoc :contig (:chr vc))
                                        (assoc :start (:start vc))
                                        (assoc :refallele (.getBaseString (:ref-allele vc)))
                                        (assoc :issubsample 0)))))))))
      (if (and subsample-params
               (> (num-variants index-file) (get-in subsample-params [:subsample :count])))
        (subsample-metrics index-file in-file ref-file subsample-params)
        (all-metrics-as-subsample index-file)))
    index-file))

(defn- maybe-fix-category
  "Provide category metrics as expected sets to match gemini usage."
  [val attr]
  (if (= :category (get-in expose-metrics [attr :x-scale :type]))
    #{val}
    val))

(defn get-raw-metrics
  "Retrieve table of raw metrics using indexed variant file"
  [in-file ref-file & {:keys [metrics use-subsample?]}]
  (let [index-db (index-variant-file in-file ref-file)
        plot-metrics (filter (partial contains? expose-metrics)
                             (or metrics (map :id (available-metrics in-file))))]
    (sql/with-connection (get-sqlite-db index-db)
      (sql/with-query-results rows
        [(str "SELECT contig, start, refallele, " (string/join ", " plot-metrics)
              " FROM metrics"
              (if use-subsample? " WHERE issubsample=1 " " "))]
              "ORDER BY contig, start"
        (doall (map (fn [orig]
                      (reduce (fn [coll x]
                                (assoc coll x (-> orig
                                                  (get (keyword (string/lower-case x)))
                                                  (maybe-fix-category x))))
                              {:id [(:contig orig) (:start orig) (:refallele orig)]}
                              plot-metrics))
                    rows))))))
