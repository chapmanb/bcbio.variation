;;Interactive functionality for scoring based web pages.

(ns bcbio.variation.score
  (:use [domina.events :only [listen! prevent-default]]
        [domina :only [log set-attr! remove-attr! swap-content!]]
        [domina.css :only [sel]])
  (:require [domina :as domina]
            [fetch.remotes :as remotes]
            [crate.core :as crate]
            [goog.dom :as dom]
            [goog.net.XhrIo :as xhr]))

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
  (swap-content! (sel (str "#" select-id))
                 (crate/html [:input {:class "input-file" :id select-id
                                      :name select-id :type "file"}])))

(defn- add-gs-input!
  "Update an input item for GenomeSpace uptake."
  [select-id placeholder]
  (swap-content! (domina/by-id select-id)
                 (crate/html
                  [:div {:class "control-group"}
                   [:label {:class "control-label" :for "gs-folder"} GenomeSpace folder]
                   [:div {:class "controls"}
                    [:select {:id "gs-folder"}
                     [:option "Test1"]
                     [:option "Test2"]]]])))

(defn ^:export upload-generalize
  "Handle generalized upload through files or GenomeSpace."
  []
  (listen! (sel "#file-choice-upload")
           :click (fn [evt]
                    (set-active-choice! "#menu-choice-upload")
                    (add-file-input! "variant-file")
                    (add-file-input! "region-file")
                    (prevent-default evt)))
  (listen! (sel "#file-choice-gs")
           :click (fn [evt]
                    (set-active-choice! "#menu-choice-gs")
                    (add-gs-input! "variant-file" "folder/filename.vcf")
                    (add-gs-input! "region-file" "folder/regions.bed")
                    (prevent-default evt))))
