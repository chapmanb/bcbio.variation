(ns bcbio.variation.workflow.xprize
  "Perform X Prize scoring workflow, handling comparison of contestant input with reference."
  (:import [java.util UUID])
  (:use [clojure.java.io]
        [bcbio.variation.compare :only [variant-comparison-from-config]]
        [bcbio.variation.combine :only [combine-variants]]
        [bcbio.variation.normalize :only [pick-best-ref]]
        [bcbio.variation.report :only [prep-scoring-table]])
  (:require [clj-yaml.core :as yaml]
            [doric.core :as doric]
            [fs.core :as fs]
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
                         (str x)
                         (:default-compare ref))]
    (->> {:dir {:out (str (fs/file (:dir work-info) "grading"))
                :prep (str (fs/file (:dir work-info) "grading" "prep"))}
          :experiments [{:sample (:sample ref)
                         :ref (:genome ref)
                         :intervals (:intervals ref)
                         :approach "grade"
                         :calls [{:name "reference"
                                  :file (:variants ref)
                                  :remove-refcalls true
                                  :normalize true}
                                 {:name "contestant"
                                  :prep true
                                  :preclean true
                                  :remove-refcalls true
                                  :ref (pick-best-ref contestant-vcf (cons (:genome ref)
                                                                           (:genome-alts ref)))
                                  :file contestant-vcf
                                  :intervals (if-let [x (:region-file work-info)]
                                               (str x)
                                               (:intervals ref))}]}]}
         yaml/generate-string
         (spit config-file))
    config-file))

(defn- write-scoring-summary
  "Output text summary file with scoring information."
  [work-info comparison]
  (let [summary-file (str (fs/file (:dir work-info)
                                   (format "%s-scoring.txt"
                                           (get-in comparison [:summary :sample]))))]
    (with-open [wtr (writer summary-file)]
      (.write wtr (str (doric/table [:metric :value]
                                    (prep-scoring-table (:metrics comparison)
                                                        (get-in comparison [:summary :sv])))
                       "\n")))
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

(defn run-scoring-analysis
  "Run scoring analysis from provided work information."
  [work-info rclient host-info config]
  (let [process-config (create-work-config work-info config)
        comparison (first (variant-comparison-from-config process-config))]
    (prepare-final-files comparison)
    {:work-info work-info
     :comparison (-> comparison
                     (assoc-in [:c-files :summary] (write-scoring-summary work-info comparison))
                     (assoc-in [:c-files :summary-html] (write-html-scoring-summary work-info comparison)))}))

(defn prep-scoring
  "Prep directory for scoring analysis."
  [params config]
  (let [tmp-dir (file (get-in config [:dir :work]) "score")
        work-id (str (UUID/randomUUID))
        cur-dir (file tmp-dir work-id)]
    (fs/mkdirs cur-dir)
    {:id work-id :dir (str cur-dir)
     :comparison-genome (:comparison-genome params)}))
