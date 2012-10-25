(ns bcbio.variation.web.process
  "Run scoring analysis, handling preparation of input files and run configuration."
  (:use [clojure.java.io]
        [ordered.map :only [ordered-map]]
        [bcbio.variation.report :only [prep-scoring-table]]
        [bcbio.variation.api.shared :only [web-config]])
  (:require [doric.core :as doric]
            [fs.core :as fs]
            [hiccup.core :as hiccup]
            [net.cgrand.enlive-html :as html]
            [ring.middleware.anti-forgery :as anti-forgery]
            [bcbio.variation.remote.core :as remote]
            [bcbio.variation.web.db :as db]))

;; ## Run scoring based on inputs from web or API

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

(defn run-scoring
  "Download form-supplied input files and start scoring analysis."
  [work-info runner rclient]
  (when-let [username (:username rclient)]
    (let [comparisons @runner]
      (db/add-analysis {:username username :files (:c-files comparisons)
                        :analysis_id (:id work-info)
                        :description (format "%s: %s" (:comparison-genome work-info)
                                             (fs/base-name (-> comparisons :exp :calls second :file)))
                        :location (:dir work-info) :type :scoring}
                       (:db @web-config)))))

(comment "Old support for upload/GenomeSpace only"
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
)
