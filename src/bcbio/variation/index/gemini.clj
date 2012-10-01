(ns bcbio.variation.index.gemini
  "Index and retrieve variant associated population genetic and disease data.
   Built on the Gemini framework: https://github.com/arq5x/gemini"
  (:use [clojure.java.io]
        [ordered.map :only [ordered-map]]
        [bcbio.variation.api.shared :only [web-config]]
        [bcbio.variation.web.db :only [get-sqlite-db]])
  (:require [clojure.java.jdbc :as sql]
            [clojure.java.shell :as shell]
            [clojure.string :as string]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]))

;; ## Variant effects

(defn- get-vep-cmd []
  (let [vep-dir (get-in @web-config [:program :vep])
        vep-file (when vep-dir (str (file (fs/expand-home vep-dir)
                                          "variant_effect_predictor.pl")))]
    (when (and vep-file (fs/exists? vep-file))
      vep-file)))

(defn run-vep
  "Run Ensembl Variant Effects Predictor on input variant file.
   Re-annotates the input file with CSQ field compatible with Gemini."
  [in-file & {:keys [re-run?]}]
  (when-let [vep-cmd (get-vep-cmd)]
    (let [out-file (itx/add-file-part in-file "vep")]
      (when (or (itx/needs-run? out-file) re-run?)
        (itx/with-tx-file [tx-out out-file]
          (shell/sh "perl" vep-cmd "-i" in-file "-o" tx-out "--vcf" "--cache"
                    "--terms" "so" "--sift" "b" "--polyphen" "b" "--hgnc" "--numbers"
                    "--fields" "Consequence,Codons,Amino_acids,Gene,HGNC,Feature,EXON,PolyPhen,SIFT")))
      out-file)))

;; ## Gemini

(defn get-gemini-cmd []
  (let [cmd (get-in @web-config [:program :gemini] "gemini")
        info (try (shell/sh cmd "-h")
                  (catch java.io.IOException _
                    {:exit -1}))]
    (when (zero? (:exit info))
      cmd)))

(defn index-variant-file
  "Pre-index a variant file with gemini"
  [in-file _ & {:keys [re-index?]}]
  (when-let [gemini-cmd (get-gemini-cmd)]
    (when in-file
      (let [index-file (str (itx/file-root in-file) "-gemini.db")]
        (when (or (itx/needs-run? index-file) re-index?)
          (let [vep-file (run-vep in-file :re-run? re-index?)]
            (itx/with-tx-file [tx-index index-file]
              (apply shell/sh
                     (concat [gemini-cmd "load" "-v"]
                             (if vep-file
                               [vep-file "-t" "VEP"]
                               [in-file])
                             [tx-index])))))
        index-file))))

(def ^{:doc "Gemini metrics to expose for query and visualization."
       :private true}
  gemini-metrics
  (ordered-map "aaf_1kg_all" {:range [0.0 1.0]
                              :desc "1000 genomes allele frequency, all populations"}
               "gms_illumina" {:range [0.0 100.0]
                               :y-scale {:type :log}
                               :desc "Genome Mappability Score with an Illumina error model"}
               "sift_score" {:range [0.0 0.1]
                             :y-scale {:type :log}
                             :desc (str "SIFT amino acid effect predictions. "
                                        "Smaller scores are more likely to be deleterious.")}
               "polyphen_score" {:range [0.5 1.0]
                                 :y-scale {:type :log}
                                 :desc (str "Polyphen amino acid effect predictions. "
                                            "Larger scores are more likely to be deleterious.")}
               "rmsk" {:viz false}))

(defn available-metrics
  "Retrieve metrics available from Gemini."
  [in-file & {:keys [include-noviz?]}]
  (let [all-metrics (->> gemini-metrics
                         (map (fn [[k v]] (assoc v :id k)))
                         (filter #(or include-noviz? (get % :viz true))))]
    (if-let [index-db (index-variant-file in-file nil)]
      (sql/with-connection (get-sqlite-db index-db)
        (letfn [(db-has-metric? [x]
                  (sql/with-query-results rows
                    [(str "SELECT chrom, start FROM variants WHERE "
                          (:id x) " IS NOT NULL LIMIT 1")]
                    (seq rows)))]
          (doall (filter db-has-metric? all-metrics))))
      all-metrics)))


(defmulti finalize-gemini-attr
  "Provide additional post-processing of gemini supplied attributes."
  (fn [attr val] (keyword (string/lower-case attr))))

(defmethod finalize-gemini-attr :sift_score
  [_ val]
  (if (nil? val) 1.0 val))

(defmethod finalize-gemini-attr :polyphen_score
  [_ val]
  (if (nil? val) 0.0 val))

(defmethod finalize-gemini-attr :gms_illumina
  [_ val]
  (if (nil? val) 100.0 val))

(defmethod finalize-gemini-attr :gms_solid
  [_ val]
  (if (nil? val) 100.0 val))

(defmethod finalize-gemini-attr :gms_iontorrent
  [_ val]
  (if (nil? val) 100.0 val))

(defmethod finalize-gemini-attr :default
  [_ val]
  val)

(defn- gemini-metric-from-row
  [row attr]
  (finalize-gemini-attr attr
                        (get row (keyword (string/lower-case attr)))))

(defn vc-attr-retriever
  "Retrieve metrics by name from a gemini index for provided VariantContexts."
  [in-file ref-file]
  (if-let [index-db (index-variant-file in-file ref-file)]
    (sql/with-connection (get-sqlite-db index-db)
      (fn [vc attr]
        (sql/with-query-results rows
          [(str "SELECT " (name attr)
                " FROM variants WHERE chrom = ? AND start = ? and ref = ?")
           (str "chr" (:chr vc)) (dec (:start vc)) (.getBaseString (:ref-allele vc))]
          (gemini-metric-from-row (first rows) attr))))
    (fn [vc attr] nil)))

(defn get-raw-metrics
  "Retrieve table of Gemini metrics keyed on variant names."
  [in-file ref-file & {:keys [metrics]}]
  (when-let [index-db (index-variant-file in-file ref-file)]
    (let [plot-metrics (filter (partial contains? gemini-metrics)
                               (or metrics (map :id (available-metrics in-file))))]
      (when (seq plot-metrics)
        (sql/with-connection (get-sqlite-db index-db)
          (sql/with-query-results rows
            [(str "SELECT chrom, start, ref, " (string/join ", " plot-metrics)
                  " FROM variants WHERE filter is NULL ORDER BY chrom, start")]
            (doall (map (fn [orig]
                          (reduce (fn [coll x]
                                    (assoc coll x (gemini-metric-from-row orig x)))
                                  {:id [(string/replace (:chrom orig) "chr" "") (inc (:start orig)) (:ref orig)]}
                                  plot-metrics))
                        rows))))))))