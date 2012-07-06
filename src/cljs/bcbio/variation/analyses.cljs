;; Interactive functionality for display of older analyses.

(ns bcbio.variation.analyses
  (:require [domina :as domina]
            [domina.css :as css]
            [domina.events :as events]
            [fetch.remotes :as remotes])
  (:require-macros [fetch.macros :as fm]))

(defn- display-selected-analysis
  "Update page with details from a selected analysis."
  [analysis-id]
  (fm/remote (get-summary analysis-id) [sum-html]
             (domina/set-html! (domina/by-id "user-analyses")
                               sum-html)))

(defn ^:export display-analyses
  "Correctly set the top level navigation toolbar."
  []
  (events/listen! (-> (domina/by-id "user-analyses")
                      (css/sel "ul")
                      domina/children)
                  :click (fn [evt]
                           (display-selected-analysis (domina/attr (events/target evt) :id))
                           (events/prevent-default evt))))
