;;Interactive functionality for scoring based web pages.

(ns bcbio.variation.score
  (:use [domina :only [log set-attr! remove-attr! swap-content!]]
        [domina.css :only [sel]])
  (:require [domina :as domina]
            [domina.events :as events]
            [fetch.remotes :as remotes]
            [crate.core :as crate]
            [goog.dom :as dom]
            [goog.net.XhrIo :as xhr])
  (:require-macros [fetch.macros :as fm]))

;; ## Display scoring results

(defn ^:export run
  "Run scoring and fetch results"
  []
  (xhr/send "/summary"
            (fn [x]
              (set! (.-innerHTML (dom/getElement "scoring-summary"))
                    (-> x .-target .getResponseText)))))

;; ## Allow multiple upload methods

(defn- set-active-choice!
  "Set the active type for our type of file upload."
  [active-id]
  (-> (sel ".nav")
      (sel ".active")
      (remove-attr! "class"))
  (-> (sel active-id)
      (set-attr! "class" "active")))

(defn- add-file-input!
  "Update an input item for file upload input."
  [select-id]
  (domina/swap-content! (domina/by-id select-id)
                        (crate/html [:input {:class "input-file" :id select-id
                                             :name select-id :type "file"}])))


(defn- update-gs-path-select!
  "Update selection box for files or directories from GenomeSpace."
  [input-id xs]
  (domina/swap-content! (domina/by-id input-id)
                        (crate/html
                         [:div]
                         [:select {:id input-id}
                          (for [x xs]
                            [:option {:value (:full x)} (:name x)])])))

(defn- update-gs-files!
  "Update file information based on parent"
  [dir file-id ftype]
  (fm/remote (list-external-files dir ftype) [files]
             (update-gs-path-select! file-id files)))

(defn- link-folders-to-files
  "Update file information as selected folder changes."
  [folder-id file-id ftype]
  (domina/log (str "Adding hook " folder-id))
  (events/listen! (domina/by-id folder-id)
                  :select (fn [evt]
                            (domina/log "Changed")
                                        ;(domina/log (events/target evt))
                            )))

(defn- add-gs-input!
  "Update an input item for GenomeSpace uptake."
  [select-id ftype]
  (let [folder-id (str "gsfolder-" select-id)
        file-id (str "gsfile-" select-id)]
    (swap-content! (domina/by-id select-id)
                   (crate/html
                    [:div {:class "control-group" :id select-id}
                     [:label {:class "control-label" :for folder-id} "Folder"]
                     [:div {:class "controls"}
                      [:input {:id folder-id :type "text"
                               :placeholder "Loading from GenomeSpace..."}]]
                     [:label {:class "control-label" :for file-id} "File"]
                     [:div {:class "controls"}
                      [:input {:id file-id :type "text"
                               :placeholder "Loading from GenomeSpace..."}]]]))
    (fm/remote (list-external-dirs) [dirs]
               (update-gs-path-select! folder-id dirs)
               (link-folders-to-files folder-id file-id ftype)
               (update-gs-files! (-> dirs first :full) file-id ftype))))

(defn ^:export upload-generalize
  "Handle generalized upload through files or GenomeSpace."
  []
  (events/listen! (sel "#file-choice-upload")
           :click (fn [evt]
                    (set-active-choice! "#menu-choice-upload")
                    (add-file-input! "variant-file")
                    (add-file-input! "region-file")
                    (events/prevent-default evt)))
  (events/listen! (sel "#file-choice-gs")
           :click (fn [evt]
                    (set-active-choice! "#menu-choice-gs")
                    (add-gs-input! "variant-file" "vcf")
                    (add-gs-input! "region-file" "bed")
                    (events/prevent-default evt))))
