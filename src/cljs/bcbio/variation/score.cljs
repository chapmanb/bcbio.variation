;;Interactive functionality for scoring based web pages.

(ns bcbio.variation.score
  (:use [domina :only [set-attr! remove-attr! swap-content!]]
        [domina.css :only [sel]])
  (:require [clojure.string :as string]
            [chosen.core :as chosen]
            [crate.core :as crate]
            [domina :as domina]
            [domina.events :as events]
            [fetch.remotes :as remotes]
            [goog.dom :as dom]
            [goog.string :as gstring]
            [goog.Timer :as timer]
            [goog.net.XhrIo :as xhr])
  (:require-macros [fetch.macros :as fm]))

;; ## Display scoring results

(defn- progress-percent
  "Rough progress points to indicate status of processing."
  [desc]
  (cond
   (gstring/startsWith desc "Starting variation") 10
   (gstring/startsWith desc "Prepare VCF, resorting to genome build: contestant") 15
   (gstring/startsWith desc "Normalize MNP and indel variants: contestant") 60
   (gstring/startsWith desc "Comparing VCFs: reference vs contestant") 75
   (gstring/startsWith desc "Summarize comparisons") 90
   (gstring/startsWith desc "Finished") 100
   :else nil))

(defn ^:export update-run-status
  "Update summary page with details about running statuses."
  [run-id]
  (fm/remote (get-status run-id) [info]
             (if (= :finished (:state info))
               (fm/remote (get-summary run-id) [sum-html]
                          (if (nil? sum-html)
                            (timer/callOnce (fn [] (update-run-status run-id)) 2000)
                             (domina/set-html! (domina/by-id "scoring-in-process")
                                               sum-html)))
               (do
                 (when-not (nil? info)
                   (domina/set-html! (domina/by-id "scoring-status")
                                     (crate/html [:p (:desc info)]))
                   (when-let [pct (progress-percent (:desc info))]
                     (domina/set-attr! (domina/by-id "scoring-progress")
                                       :style (str "width: " pct "%"))))
                 (timer/callOnce (fn [] (update-run-status run-id)) 2000)))))

;; ## Allow multiple upload methods

(defn- set-active-choice!
  "Set the active type for our type of file upload."
  [active-id]
  (-> (domina/by-id "upload-choices")
      (sel ".active")
      (remove-attr! :class))
  (-> (sel active-id)
      (set-attr! :class "active")))

(defn- add-file-input!
  "Update an input item for file upload input."
  [select-id]
  (domina/swap-content! (domina/by-id select-id)
                        (crate/html [:input {:class "input-file" :id select-id
                                             :name select-id :type "file"}])))

(defn- gs-paths-to-chosen [xs]
  (map (fn [x] {:value (:full x) :text (:name x)}) xs))

(defn- update-gs-files!
  "Update file information based on parent"
  [file-chosen file-id dir ftype]
  (let [final-form-id (string/join "-" (cons "gs" (rest (string/split file-id #"-"))))]
    (fm/remote (list-external-files dir ftype) [files]
               (chosen/options file-chosen (gs-paths-to-chosen files))
               (domina/set-value! (domina/by-id final-form-id) (chosen/selected file-chosen))
               (add-watch file-chosen :change
                          (fn [fname]
                            (domina/set-value! (domina/by-id final-form-id) fname))))))

(defn- add-gs-input!
  "Update an input item for GenomeSpace uptake."
  [select-id ftype]
  (let [folder-id (str "gsfolder-" select-id)
        file-id (str "gsfile-" select-id)]
    (swap-content! (domina/by-id select-id)
                   (crate/html
                    [:div {:id select-id}
                     [:select {:id folder-id :data-placeholder "GenomeSpace Folder"}]
                     [:select {:id file-id :data-placeholder "GenomeSpace File"}]]))
    (let [folder-chosen (chosen/ichooseu! (str "#" folder-id))
          file-chosen (chosen/ichooseu! (str "#" file-id))]
      (fm/remote (list-external-dirs) [dirs]
                 (chosen/options folder-chosen (gs-paths-to-chosen dirs))
                 (when-let [cur-dir (chosen/selected folder-chosen)]
                   (update-gs-files! file-chosen file-id cur-dir ftype))
                 (add-watch folder-chosen :change
                            (fn [dir]
                              (update-gs-files! file-chosen file-id dir ftype)))))))

(defn- prep-genome-selector
  "Prepare genome selector to pick analysis genome."
  []
  (let [genome-chosen (chosen/ichooseu! "#comparison-genome")]
    (fm/remote (get-genomes) [genomes]
               (chosen/options genome-chosen genomes))))

(defn- set-upload-active []
  (set-active-choice! "#menu-choice-upload")
  (add-file-input! "variant-file")
  (add-file-input! "region-file"))

(defn- set-gs-active []
  (set-active-choice! "#menu-choice-gs")
  (add-gs-input! "variant-file" "vcf")
  (add-gs-input! "region-file" "bed"))

(defn ^:export set-navigation
  "Correctly set the active top level navigation toolbar."
  []
  (let [loc (-> (.toString window.location ())
                (string/split #"/")
                last)]
    (doseq [list-item (domina/children (domina/by-id "top-navbar"))]
      (if (= (str "/" loc)
             (-> (domina/children list-item)
                 first
                 (domina/attr :href)))
        (domina/set-attr! list-item :class "active")
        (domina/remove-attr! list-item :class)))))

(defn ^:export upload-generalize
  "Handle generalized upload through files or GenomeSpace."
  []
  (prep-genome-selector)
  (fm/remote (get-username) [user]
             (when-not (nil? user)
               (set-gs-active)))
  (events/listen! (sel "#file-choice-upload")
                  :click (fn [evt]
                           (set-upload-active)
                           (events/prevent-default evt)))
  (events/listen! (sel "#file-choice-gs")
                  :click (fn [evt]
                           (set-gs-active)
                           (events/prevent-default evt))))
