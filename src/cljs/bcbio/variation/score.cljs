;;Interactive functionality for scoring based web pages.

(ns bcbio.variation.score
  (:require [goog.dom :as dom]
            [goog.net.XhrIo :as xhr]))

(defn ^:export run
  "Run scoring and fetch results"
  []
  (xhr/send "/summary"
            (fn [x]
              (set! (.-innerHTML (dom/getElement "scoring-summary"))
                    (-> x .-target .getResponseText)))))
