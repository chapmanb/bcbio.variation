;; Handle login and logout management.

(ns bcbio.variation.login
  (:use [domina.events :only [listen! prevent-default]]
        [domina.css :only [sel]]
        [domina.xpath :only [xpath]])
  (:require [domina :as domina]
            [fetch.remotes :as remotes]
            [crate.core :as crate])
  (:require-macros [fetch.macros :as fm]))

(defn- form-input [content]
  [(keyword (domina/attr content "name"))
   (domina/value content)])

(defn- logged-in-html [user]
  (crate/html
   [:div {:class "btn-group"}
    [:button {:class "btn btn-info dropdown-toggle" :data-toggle "dropdown"}
     [:i {:class "icon-user icon-white" :style "margin-right: 6px"}]
     user
     [:span {:class "caret" :style "margin-left: 6px"}]]
    [:ul {:class "dropdown-menu"}
     [:li [:a {:id "logout-btn" :href "#"} "Logout"]]]]))

(defn- logged-out-html []
  (crate/html
   [:form {:class "navbar-form"}
    [:input {:type "text" :name "username" :class "input-medium" :placeholder "GenomeSpace Username"
             :style "margin-top: 2px; margin-right: 4px"}]
    [:input {:type "password" :name "password" :class "input-medium" :placeholder "Password"
             :style "margin-top: 2px;"}]
    [:button {:type "submit" :id "login-btn" :class "btn-info btn"} "Login"]
    [:a {:class "btn btn-success" :href "http://www.genomespace.org/register" :target "_blank"}
     "Register"]]))

(declare login-listeners)

(defn- update-login
  "Check for logged in users, updating user management region accordingly."
  []
  (fm/remote (get-username) [user]
             (domina/set-html! (domina/by-id "user-manage")
                   (if (nil? user)
                     (logged-out-html)
                     (logged-in-html user)))
             (login-listeners)))

(defn- login-listeners
  "Add listeners for login related clicks. Updated when DOM changes."
  []
  (listen! (domina/by-id "login-btn")
           :click (fn [evt]
                    (let [login-vals (->> (xpath "//div[@id='user-manage']/form/input")
                                          (domina/nodes)
                                          (map form-input)
                                          (into {}))]
                      (fm/remote (login login-vals) [result]
                                 (if (nil? result)
                                   (js/alert "Invalid username/password")
                                   (update-login))))
                    (prevent-default evt)))
  (listen! (domina/by-id "logout-btn")
           :click (fn [evt]
                    (fm/remote (logout) []
                               (update-login))
                    (prevent-default evt))))

(defn ^:export handle-login
  "Catch login details, check login and update header."
  []
  (update-login))
