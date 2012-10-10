(ns bcbio.variation.web.process
  "Run scoring analysis, handling preparation of input files and run configuration."
  (:import [java.util.UUID])
  (:use [clojure.java.io]
        [ordered.map :only [ordered-map]]
        [bcbio.variation.combine :only [combine-variants]]
        [bcbio.variation.compare :only [variant-comparison-from-config]]
        [bcbio.variation.normalize :only [pick-best-ref]]
        [bcbio.variation.report :only [prep-scoring-table]]
        [bcbio.variation.api.shared :only [web-config]])
  (:require [clojure.string :as string]
            [clj-yaml.core :as yaml]
            [doric.core :as doric]
            [fs.core :as fs]
            [hiccup.core :as hiccup]
            [net.cgrand.enlive-html :as html]
            [ring.middleware.anti-forgery :as anti-forgery]
            [bcbio.variation.remote.core :as remote]
            [bcbio.variation.web.db :as db]))

;; ## Run scoring based on inputs from web or API

(defn create-work-config
  "Create configuration for processing inputs using references supplied in config."
  [work-info config]
  (if-not (fs/exists? (:dir work-info))
    (fs/mkdirs (:dir work-info)))
  (let [config-file (str (fs/file (:dir work-info) "process.yaml"))
        ref (first (filter #(= (:sample %) (:comparison-genome work-info))
                           (:ref config)))
        contestant-vcf (if-let [x (get-in work-info [:in-files :variant-file])]
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
                                  :intervals (if-let [x (get-in work-info [:in-files :region-file])]
                                               (str x)
                                               (:intervals ref))}]}]}
         yaml/generate-string
         (spit config-file))
    config-file))

(defn html-summary-table
  "Generate a summary table of scoring results."
  [comparisons]
  (let [scoring-table (prep-scoring-table (:metrics comparisons)
                                          (get-in comparisons [:summary :sv]))]
    (apply str
           (-> (str (doric/table ^{:format doric/html} [:metric :value] scoring-table))
               java.io.StringReader.
               html/html-resource
               (html/transform [:table] (html/set-attr :class "table table-condensed"))
               html/emit*))))

(defn- base-page-w-content
  "Update main page HTML with specified content at selector"
  [selector new-hiccup-html]
  (-> (html/html-resource (file (get-in @web-config [:dir :html-root] "public") "index.html"))
      (html/transform selector (-> new-hiccup-html
                                   java.io.StringReader.
                                   html/html-resource
                                   html/content))
      html/emit*
      (#(apply str %))))

(defn html-submit-page
  "Return main HTML page with submit information."
  []
  (base-page-w-content [:div#anti-forgery]
                       (hiccup/html [:input {:class "hidden"
                                             :name "__anti-forgery-token"
                                             :value anti-forgery/*anti-forgery-token*}])))

(defn html-scoring-summary
  "Generate summary of scoring results for display."
  [comparisons run-id]
  (hiccup/html
   [:h3 "Summary"]
   [:div {:id "score-table"}
    (html-summary-table comparisons)]
   [:h3 "Variant files in VCF format"]
   [:div {:id "variant-file-download"}
    [:ul
     (for [[key txt] [["concordant" "Concordant variants"]
                      ["discordant" "Discordant variants"]
                      ["discordant-missing" "Missing variants"]
                      ["phasing" "Variants with phasing errors"]]]
       [:li [:a {:href (format "/scorefile/%s/%s" run-id key)} txt]])]]))

(defn scoring-html
  "Update main page HTML with content for scoring."
  [run-id]
  (let [html-dir (get-in @web-config [:dir :html-root] "public")
        template-dir (str (fs/file html-dir "template"))]
    (base-page-w-content
     [:div#main-content]
     (hiccup/html
      [:div {:id "scoring-in-process"}
       [:h3 "Status"]
       [:div {:id "scoring-status"} "Downloading input files"]
       [:div {:class "progress"}
        [:div {:id "scoring-progress"
               :class "bar" :style "width: 0%"}]]
       (slurp (fs/file template-dir "scoring.html"))
       [:script {:src "js/score.js"}]
       [:script (format "bcbio.variation.score.update_run_status('%s');"
                        run-id)]]))))

(defn analyses-html
  "Update main page with list of performed analyses."
  [username]
  (base-page-w-content
   [:div#main-content]
   (hiccup/html
    (if (nil? username)
      [:p "Please login to display previously run analyses."]
      [:div {:id "user-analyses" :class "container"}
       [:h3 "Previous analyses"]
       [:ul {:class "nav nav-tabs nav-stacked"}
        (map (fn [x]
               [:li
                [:a {:href "#" :id (:analysis_id x)}
                 (format "%s -- %s" (:description x)
                         (-> (java.text.SimpleDateFormat. "dd MMM yyyy HH:mm" )
                             (.format (:created x))))]])
             (db/get-analyses username :scoring (:db @web-config)))]
       [:script {:src "js/score.js"}]
       [:script "bcbio.variation.analyses.display_analyses()"]]))))

(defn upload-results
  "Upload output files to GenomeSpace."
  [rclient work-info comparison]
  (let [out-files (map #(get-in comparison [:c-files %])
                       [:concordant :discordant :discordant-missing :phasing-error])
        summary-file (str (fs/file (:dir work-info)
                                   (format "%s-scoring.txt"
                                           (get-in comparison [:summary :sample]))))]
    (with-open [wtr (writer summary-file)]
      (.write wtr (str (doric/table [:metric :value]
                                    (prep-scoring-table (:metrics comparison)
                                                        (get-in comparison [:summary :sv])))
                       "\n")))
    (let [subdir "xprize"
          gs-dir (str (fs/file (fs/parent (last (string/split (:gs-variant-file work-info) #":" 2))) subdir))]
      (doseq [x (cons summary-file out-files)]
        (remote/put-file rclient gs-dir x {})))))

(defn- prepare-final-files
  "Merge standard and structural variant outputs into final set of upload files."
  [comparisons]
  (letfn [(merge-files-into [comparisons orig-kw addin-kw]
            (let [ref (get-in comparisons [:exp :ref])
                  orig (get-in comparisons [:c-files orig-kw])
                  addin (get-in comparisons [:c-files addin-kw])
                  combine-vcf (combine-variants [orig addin] ref :merge-type :full
                                                :quiet-out? true)]
              (fs/rename combine-vcf orig)
              (fs/rename (str combine-vcf ".idx") (str orig ".idx"))))]
      (merge-files-into comparisons :concordant :sv-concordant)
      (merge-files-into comparisons :discordant :sv-contestant-discordant)
      (merge-files-into comparisons :discordant-missing :sv-reference-discordant)))

(defn run-scoring-analysis
  "Run scoring analysis from provided work information."
  [work-info rclient]
  (let [process-config (create-work-config work-info @web-config)
        comparisons (first (variant-comparison-from-config process-config))]
    (prepare-final-files comparisons)
    (when-not (or (nil? (:conn rclient)) (nil? (:gs-variant-file work-info)))
      (upload-results rclient work-info comparisons))
    (spit (file (:dir work-info) "scoring-summary.html")
          (html-scoring-summary comparisons (:id work-info)))
    comparisons))

(defmulti get-input-files
  "Prepare working directory, downloading input files."
  (fn [work-info params _]
    (let [gs-info (:gs-variant-file params)]
      (if (or (nil? gs-info) (empty? gs-info)) :upload :gs))))

(defmethod get-input-files :upload
  [work-info params _]
  (letfn [(download-file [tmp-dir params kw]
            (let [cur-param (get params kw)
                  out-file (fs/file tmp-dir (:filename cur-param))]
              [kw (when (> (:size cur-param) 0)
                    (copy (:tempfile cur-param) out-file)
                    (str out-file))]))]
    (into {} (map (partial download-file (:dir work-info) params)
                  [:variant-file :region-file]))))

(defmethod get-input-files :gs
  [work-info params rclient]
  (letfn [(gs-download [tmp-dir params kw]
            (let [gs-file (get params (keyword (str "gs-" (name kw))))]
              [kw (when-not (or (nil? gs-file) (empty? gs-file))
                    (remote/get-file gs-file rclient :out-dir tmp-dir))]))]
    (into {} (map (partial gs-download (:dir work-info) params)
                  [:variant-file :region-file]))))

(defn run-scoring
  "Download form-supplied input files and start scoring analysis."
  [orig-work-info params rclient]
  (letfn [(add-gs-info [work-info params]
            (let [remote-fname (:gs-variant-file params)]
              (if-not (or (nil? remote-fname) (empty? remote-fname))
                (assoc work-info :gs-variant-file remote-fname)
                work-info)))]
    (let [work-info (-> orig-work-info
                        (add-gs-info params)
                        (assoc :in-files (get-input-files orig-work-info params
                                                          rclient)))
          comparisons (run-scoring-analysis work-info rclient)]
      (when-let [username (:username rclient)]
        (db/add-analysis {:username username :files (:c-files comparisons)
                          :analysis_id (:id work-info)
                          :description (format "%s: %s" (:comparison-genome work-info)
                                               (fs/base-name (-> comparisons :exp :calls second :file)))
                          :location (:dir work-info) :type :scoring}
                         (:db @web-config))))))

(defn prep-scoring
  "Prep directory for scoring analysis."
  [params]
  (letfn [(prep-tmp-dir []
            (let [tmp-dir (file (get-in @web-config [:dir :work]) "score")
                  work-id (str (java.util.UUID/randomUUID))
                  cur-dir (file tmp-dir work-id)]
              (fs/mkdirs cur-dir)
              {:id work-id :dir (str cur-dir)
               :comparison-genome (:comparison-genome params)}))]
    (let [work-info (prep-tmp-dir)]
      {:work-info work-info
       :out-html (scoring-html (:id work-info))})))

;; ## File retrieval from processing

(defn- get-run-info
  "Retrieve run information from stored database or current work-info."
  [run-id username session-work-info]
  (if (nil? username)
    [(:comparison-genome session-work-info) (:dir session-work-info)]
    (let [work-info (->> (db/get-analyses username :scoring (:db @web-config))
                         (filter #(= run-id (:analysis_id %)))
                         first)
          sample-name (when-not (nil? work-info)
                        (first (string/split (:description work-info) #":")))]
      [sample-name (:location work-info)])))

(defn get-variant-file
  "Retrieve processed output file for web display."
  [run-id name username work-info]
  (letfn [(sample-file [sample-name ext]
            (let [base-name "contestant-reference"]
              (format "%s-%s-%s" sample-name base-name ext)))]
    (let [[sample-name base-dir] (get-run-info run-id username work-info)
          file-map {"concordant" (sample-file sample-name "concordant.vcf")
                    "discordant" (sample-file sample-name "discordant.vcf")
                    "discordant-missing" (sample-file sample-name "discordant-missing.vcf")
                    "phasing" (sample-file sample-name "phasing-error.vcf")}
          work-dir (when-not (nil? base-dir) (fs/file base-dir "grading"))
          name (get file-map name)
          fname (if-not (or (nil? work-dir)
                            (nil? name)) (str (fs/file work-dir name)))]

      (if (and (not (nil? fname)) (fs/exists? fname))
        (slurp fname)
        "Variant file not found"))))
