(ns bcbio.variation.workflow.xprize
  "Perform X Prize scoring workflow, handling comparison of contestant input with reference."
  (:import [java.util UUID])
  (:use [clojure.java.io]
        [bcbio.variation.config :only [traceback-to-log load-config]]
        [bcbio.variation.compare :only [variant-comparison-from-config]]
        [bcbio.variation.combine :only [combine-variants]]
        [bcbio.variation.normalize :only [pick-best-ref]]
        [bcbio.variation.report :only [prep-scoring-table]])
  (:require [clj-yaml.core :as yaml]
            [doric.core :as doric]
            [me.raynes.fs :as fs]
            [hiccup.core :as hiccup]
            [net.cgrand.enlive-html :as html]))

(defn create-work-config
  "Create configuration for processing inputs using references supplied in config."
  [work-info config]
  (if-not (fs/exists? (:dir work-info))
    (fs/mkdirs (:dir work-info)))
  (let [config-file (str (fs/file (:dir work-info) "process.yaml"))
        ref (first (filter #(= (:sample %) (:comparison-genome work-info))
                           (:ref config)))
        contestant-vcf (if-let [x (:variant-file work-info)]
                         (str (fs/file x))
                         (:default-compare ref))]
    (->> {:dir {:out (str (fs/file (:dir work-info) "grading"))
                :prep (str (fs/file (:dir work-info) "grading" "prep"))}
          :experiments [{:sample (:sample ref)
                         :ref (:genome ref)
                         :intervals (:intervals ref)
                         :approach "grade"
                         :calls [{:name "reference"
                                  :file (:variants ref)
                                  :normalize false}
                                 {:name "contestant"
                                  :prep true
                                  :preclean true
                                  :remove-refcalls true
                                  :ref (pick-best-ref contestant-vcf (cons (:genome ref)
                                                                           (:genome-alts ref)))
                                  :file contestant-vcf
                                  :intervals (if-let [x (:region-file work-info)]
                                               (str (fs/file x))
                                               (:intervals ref))}]}]}
         yaml/generate-string
         (spit config-file))
    config-file))

(defn- write-scoring-summary
  "Output text summary file with scoring information."
  [work-info comparison]
  (let [summary-file (str (fs/file (:dir work-info)
                                   (format "%s-scoring.tsv"
                                           (get-in comparison [:summary :sample]))))]
    (with-open [wtr (writer summary-file)]
      (doseq [x (prep-scoring-table (:metrics comparison)
                                    (get-in comparison [:summary :sv]))]
        (.write wtr (format "%s\t%s\n" (:metric x) (:value x)))))
    summary-file))

(defn- html-summary-table
  "Generate a summary table of scoring results."
  [comparison]
  (let [scoring-table (prep-scoring-table (:metrics comparison)
                                          (get-in comparison [:summary :sv]))]
        
    (apply str
           (-> (str (doric/table ^{:format doric/html} [:metric :value] scoring-table))
               java.io.StringReader.
               html/html-resource
               (html/transform [:table] (html/set-attr :class "table table-condensed"))
               html/emit*))))

(defn- write-html-scoring-summary
  "Generate summary of scoring results for display."
  [work-info comparison]
  (let [out-file (str (file (:dir work-info) "scoring-summary.html"))]
    (spit out-file
          (hiccup/html
           [:h3 "Summary"]
           [:div {:id "score-table"}
            (html-summary-table comparison)]
           [:h3 "Variant files in VCF format"]
           [:div {:id "variant-file-download"}
            [:ul
             (for [[key txt] [["concordant" "Concordant variants"]
                              ["discordant" "Discordant variants"]
                              ["discordant-missing" "Missing variants"]
                              ["phasing" "Variants with phasing errors"]]]
               [:li [:a {:href (format "/dataset/%s/%s" (:id work-info) key)} txt]])]]))
    out-file))

(defn- prepare-final-files
  "Merge standard and structural variant outputs into final set of upload files."
  [comparison]
  (letfn [(merge-files-into [comparison orig-kw addin-kw]
            (let [ref (get-in comparison [:exp :ref])
                  orig (get-in comparison [:c-files orig-kw])
                  addin (get-in comparison [:c-files addin-kw])
                  combine-vcf (combine-variants [orig addin] ref :merge-type :full
                                                :quiet-out? true)]
              (fs/rename combine-vcf orig)
              (fs/rename (str combine-vcf ".idx") (str orig ".idx"))))]
      (merge-files-into comparison :concordant :sv-concordant)
      (merge-files-into comparison :discordant :sv-contestant-discordant)
      (merge-files-into comparison :discordant-missing :sv-reference-discordant)))

(defn run-scoring-analysis*
  "Run X Prize scoring analysis from provided work information."
  [work-info rclient config-file]
  (let [comparison (first (variant-comparison-from-config config-file))]
    (prepare-final-files comparison)
    {:work-info work-info
     :comparison (-> comparison
                     (assoc-in [:c-files :summary] (write-scoring-summary work-info comparison))
                     (assoc-in [:c-files :summary-html] (write-html-scoring-summary work-info comparison)))}))

(defn run-scoring-analysis
  "Safe running of X Prize workflow with exception catching."
  [work-info rclient input-config]
  (prn input-config)
  (let [config-file (create-work-config work-info input-config)]
    (try
      (run-scoring-analysis* work-info rclient config-file)
      (catch Exception e
        (do
          (traceback-to-log e (load-config config-file))
          (throw e))))))

(defn prep-scoring
  "Prep directory for scoring analysis."
  [params config]
  (let [tmp-dir (file (get-in config [:dir :work]) "score")
        work-id (str (UUID/randomUUID))
        cur-dir (file tmp-dir work-id)]
    (fs/mkdirs cur-dir)
    (-> params
        (assoc :id work-id)
        (assoc :dir (str cur-dir)))))
