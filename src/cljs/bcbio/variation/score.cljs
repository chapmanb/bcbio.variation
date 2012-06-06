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
