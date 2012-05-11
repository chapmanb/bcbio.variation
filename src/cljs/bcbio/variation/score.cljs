;;Interactive functionality for scoring based web pages.

(ns bcbio.variation.score
  (:use [domina.events :only [listen! prevent-default]]
        [domina :only [log set-attr! remove-attr! swap-content!]]
        [domina.css :only [sel]])
  (:require [goog.dom :as dom]
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
                 (str "<input class='input-file' id='" select-id "' "
                      "name='" select-id "' type='file'></input>")))

(defn- add-s3-input!
  "Update an input item for S3 uptake."
  [select-id placeholder]
  (-> (sel (str "#" select-id))
      (set-attr! "type" "text")
      (set-attr! "placeholder" placeholder)
      (remove-attr! "class")))

(defn ^:export upload-generalize
  "Handle generalized upload through files or S3."
  []
  (listen! (sel "#file-choice-upload")
           :click (fn [evt]
                    (set-active-choice! "#menu-choice-upload")
                    (add-file-input! "variant-file")
                    (add-file-input! "region-file")
                    (prevent-default evt)))
  (listen! (sel "#file-choice-s3")
           :click (fn [evt]
                    (set-active-choice! "#menu-choice-s3")
                    (add-s3-input! "variant-file" "s3://bucketname/filename.vcf")
                    (add-s3-input! "region-file" "s3://bucketname/regions.bed")
                    (prevent-default evt))))
