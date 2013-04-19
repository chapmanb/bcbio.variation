(ns bcbio.variation.index.gemini
  "Index and retrieve variant associated population genetic and disease data.
   Built on the Gemini framework: https://github.com/arq5x/gemini"
  (:use [clojure.java.io]
        [ordered.map :only [ordered-map]]
        [bcbio.variation.api.shared :only [web-config]]
        [bcbio.variation.web.db :only [get-sqlite-db get-sqlite-db-pool]])
  (:require [clojure.java.jdbc :as sql]
            [clojure.java.shell :as shell]
            [clojure.string :as string]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]))

;; ## Gemini

(defn get-gemini-cmd []
  (let [cmd (get-in @web-config [:program :gemini] "gemini")
        info (try (shell/sh cmd "-h")
                  (catch java.io.IOException _
                    {:exit -1}))]
    (when (zero? (:exit info))
      cmd)))

(defn- has-snpeff-anns?
  "Check if the input file contains snpEff annotations."
  [in-file]
  (with-open [rdr (reader in-file)]
    (->> (line-seq rdr)
         (take-while #(.startsWith % "##"))
         (filter #(.startsWith % "##SnpEff"))
         count
         pos?)))

(defn index-variant-file
  "Pre-index a variant file with gemini, handling snpEff annotations."
  [in-file _ & {:keys [re-index?]}]
  (when-let [gemini-cmd (get-gemini-cmd)]
    (when in-file
      (let [index-file (str (itx/file-root in-file) "-gemini.db")]
        (when (or (itx/needs-run? index-file) re-index?)
          (itx/with-tx-file [tx-index index-file]
            (apply shell/sh
                   (concat [gemini-cmd "load" "-v" in-file]
                           (when (has-snpeff-anns? in-file)
                             ["-t" "snpEff"])
                           [tx-index]))))
        index-file))))

(def ^{:doc "Gemini metrics to expose for query and visualization."
       :private true}
  gemini-metrics
  (ordered-map
   "aaf_1kg_all" {:range [0.0 1.0]
                  :desc "1000 genomes allele frequency, all populations"}
   "gms_illumina" {:range [0.0 100.0]
                   :y-scale {:type :log}
                   :desc "Genome Mappability Score with an Illumina error model"}
   "in_cse" {:x-scale {:type :category}
             :desc "Presence of variant in an error prone genomic position"}
   "rmsk" {:x-scale {:type :category}
           :desc "Repeat status: is the variant in a known repeat region"}
   "type" {:x-scale {:type :category}
           :rows {:type "" :sub_type ""}
           :desc "Type of variant change"}
   "zygosity" {:x-scale {:type :category}
               :rows {:num_hom_ref "homozygous ref"
                      :num_het "heterozygous"
                      :num_hom_alt "homozygous"}
               :desc "Allele types present in individuals"}
   "encode_consensus_gm12878" {:x-scale {:type :category}
                               :desc "Chromatin status: consensus from ENCODE for NA12878"}
   "in_public" {:x-scale {:type :category}
                :rows {:in_dbsnp "dbSNP"
                       :in_hm3 "HapMap3"
                       :in_esp "ESP"
                       :in_1kg "1000genomes"}
                :desc "Presence in large variant projects like dbSNP and 1000 genomes"}
   "is_coding" {:x-scale {:type :category}
                :desc "Type of coding transcript influenced by variant"}
   "impact_severity" {:x-scale {:type :category}
                      :desc "Severity of variant impact on coding region"}))

(defn- attr->colnames
  "Convert a gemini attribute into potentially multiple gemini column names."
  [attr]
  (->> (if-let [rows (get-in gemini-metrics [(name attr) :rows])]
         (keys rows)
         [attr])
       (map name)))

(defn available-metrics
  "Retrieve metrics available from Gemini."
  [in-file & {:keys [noviz?]}]
  (let [all-metrics (->> gemini-metrics
                         (map (fn [[k v]] (assoc v :id k)))
                         (filter #(or noviz? (get % :viz true))))]
    (if-let [index-db (index-variant-file in-file nil)]
      (sql/with-connection (get-sqlite-db index-db)
        (letfn [(db-has-metric? [x]
                  (sql/with-query-results rows
                    [(str "SELECT chrom, start FROM variants WHERE "
                          (->> (attr->colnames (:id x))
                               (map #(str % " IS NOT NULL"))
                               (string/join " OR "))
                           " LIMIT 1")]
                    (seq rows)))]
          (doall (filter db-has-metric? all-metrics))))
      all-metrics)))


(defmulti finalize-gemini-attr
  "Provide additional post-processing of gemini supplied attributes."
  (fn [attr row] (keyword (string/lower-case attr))))

(defmethod finalize-gemini-attr :sift_score
  [_ row]
  (let [val (first (vals row))]
    (if (nil? val) 1.0 val)))

(defmethod finalize-gemini-attr :polyphen_score
  [_ row]
  (let [val (first (vals row))]
    (if (nil? val) 0.0 val)))

(defmethod finalize-gemini-attr :gms_illumina
  [_ row]
  (let [val (first (vals row))]
    (if (nil? val) 100.0 val)))

(defmethod finalize-gemini-attr :gms_solid
  [_ row]
  (let [val (first (vals row))]
    (if (nil? val) 100.0 val)))

(defmethod finalize-gemini-attr :gms_iontorrent
  [_ row]
  (let [val (first (vals row))]
    (if (nil? val) 100.0 val)))

(defmethod finalize-gemini-attr :in_cse
  [_ row]
  (let [val (first (vals row))]
    (if (and (not (nil? val)) (pos? val)) #{"error-prone"} #{"standard"})))

(defmethod finalize-gemini-attr :rmsk
  [_ row]
  (let [val (first (vals row))]
    #{(if (nil? val) "non-repeat" "repeat")}))

(defmethod finalize-gemini-attr :type
  [_ row]
  (set (map #(case %
               "ts" "transition"
               "tv" "transversion"
               "ins" "insertion"
               "del" "deletion"
               %)
            (vals row))))

(defn- row->names
  "Convert a row into pre-configured names based on gemini-metrics :rows"
  [attr row]
  (let [row-names (get-in gemini-metrics [attr :rows])]
    (reduce (fn [coll [k v]]
              (if (and (not (nil? v))
                       (pos? v))
                (conj coll (get row-names (keyword k)))
                coll))
            #{} row)))

(defmethod finalize-gemini-attr :zygosity
  [attr row]
  (row->names attr row))

(defmethod finalize-gemini-attr :encode_consensus_gm12878
  ^{:doc "ENCODE chromatin segment predictions, from Table 3 of doi:10.1038/nature11247"}
  [_ row]
  (let [val (first (vals row))]
    #{(case val
        "CTCF" "CTCF-enriched"
        "E" "Enhancer"
        "PF" "Promoter flanking"
        "R" "Repressed"
        "TSS" "Promoter with TSS"
        "T" "Transcribed"
        "WE" "Weak enchancer"
        "Unknown")}))

(defmethod finalize-gemini-attr :in_public
  [attr row]
  (let [publics (row->names attr row)]
    (if (seq publics)
      publics
      #{"unique"})))

(defmethod finalize-gemini-attr :is_coding
  [_ row]
  (let [val (first (vals row))]
    (if (and (not (nil? val)) (pos? val)) #{"coding"} #{"noncoding"})))

(defmethod finalize-gemini-attr :impact_severity
  [_ row]
  (let [val (first (vals row))]
    #{val}))

(defmethod finalize-gemini-attr :default
  [_ row]
  (let [val (first (vals row))]
    val))

(defn- gemini-metric-from-row
  [row attr]
  (let [colnames (attr->colnames attr)]
    (finalize-gemini-attr attr
                          (zipmap colnames
                                  (map #(get row (keyword (string/lower-case %))) colnames)))))


(defn vc-attr-retriever
  "Retrieve metrics by name from a gemini index for provided VariantContexts."
  [in-file ref-file]
  (if-let [index-db (index-variant-file in-file ref-file)]
    (let [pool (get-sqlite-db-pool index-db)]
      (fn [vc attr]
        (sql/with-connection pool
          (sql/with-query-results rows
            [(str "SELECT " (string/join "," (attr->colnames attr))
                  " FROM variants WHERE chrom = ? AND start = ? and ref = ?")
             (str "chr" (:chr vc)) (dec (:start vc)) (.getBaseString (:ref-allele vc))]
            (gemini-metric-from-row (first rows) attr)))))
    (fn [vc attr] nil)))

(defn get-raw-metrics
  "Retrieve table of Gemini metrics keyed on variant names."
  [in-file ref-file & {:keys [metrics noviz?]}]
  (when-let [index-db (index-variant-file in-file ref-file)]
    (let [plot-metrics (filter (partial contains? gemini-metrics)
                               (or metrics (map :id (available-metrics in-file
                                                                       :noviz? noviz?))))]
      (when (seq plot-metrics)
        (sql/with-connection (get-sqlite-db index-db)
          (sql/with-query-results rows
            [(str "SELECT chrom, start, ref, "
                  (string/join ", " (mapcat attr->colnames plot-metrics))
                  " FROM variants WHERE (filter is NULL or filter = 'PASS')"
                  " ORDER BY chrom, start")]
            (doall (map (fn [orig]
                          (reduce (fn [coll x]
                                    (assoc coll x (gemini-metric-from-row orig x)))
                                  {:id [(string/replace (:chrom orig) "chr" "") (inc (:start orig)) (:ref orig)]}
                                  plot-metrics))
                        rows))))))))
