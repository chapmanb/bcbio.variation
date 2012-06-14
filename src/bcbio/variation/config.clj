(ns bcbio.variation.config
  "Load and prepare inputs from YAML configuration files."
  (:use [clojure.java.io])
  (:require [clj-logging-config.log4j :as log4j]
            [clojure.tools.logging :as log]
            [clojure.tools.logging.impl :as log-impl]
            [clj-yaml.core :as yaml]
            [fs.core :as fs]
            [pallet.algo.fsm.fsm :as fsm-base]
            [pallet.algo.fsm.fsm-dsl :as fsm]
            [pallet.algo.fsm.event-machine :as event-machine]))

;; ## Logging

(defn prep-comparison-fsm
  "Define a finite state machine of transitions during comparison processes."
  []
  (letfn [(log-transition [_ new-state]
            (let [out (format "State %s => %s" (:state-kw new-state) (:state-data new-state))]
              (println out)
              (log/log :info out)))]
    (log-transition nil {:state-kw :begin :state-data {:desc "Starting variation analysis"}})
    (event-machine/event-machine
     (fsm/event-machine-config
      (fsm/using-fsm-features (fsm-base/with-transition-observer log-transition))
      (fsm/initial-state :begin)
      (fsm/initial-state-data {})
      (fsm/state :begin
                 (fsm/valid-transitions :merge :clean))
      (fsm/state :clean
                 (fsm/valid-transitions :merge :prep))
      (fsm/state :prep
                 (fsm/valid-transitions :normalize))
      (fsm/state :normalize
                 (fsm/valid-transitions :merge))
      (fsm/state :merge)))))

(defn configure-logging
  "Setup output file logging based on configuration"
  [config]
  (alter-var-root (var log/*logger-factory*) (constantly (log-impl/log4j-factory)))
  (let [out-file (file (get-in config [:dir :out]) "processing-status.log")
        pattern "%d [%-5p] %m%n"]
    (log4j/set-loggers! :config {:level :info}
                        :root {:level :info :pattern pattern :out :console}
                        "bcbio.variation.config" {:out out-file :level :info :pattern pattern}
                        "org.broadinstitute.sting.gatk.refdata.tracks.RMDTrackBuilder"
                        {:level :warn}))
  (let [fsm (prep-comparison-fsm)]
    (-> config
        (assoc :fsm fsm)
        (assoc :transition (fn [state desc]
                             ((:transition fsm) #(assoc % :state-kw state
                                                        :state-data {:desc desc})))))))

;; ## Configuration

(defn- add-dir-files
  "Add files of interest in a directory with the given extension.
  This allows batch processing of directories."
  [config exts]
  (letfn [(files-from-dir [dir]
            (->> (fs/list-dir dir)
                 (filter #(contains? exts (fs/extension %)))
                 (map #(str (fs/file dir %)))))
          (process-call [call]
            (if-let [dir (:dir call)]
              (assoc call :file (files-from-dir dir))
              call))
          (process-exp [exp]
            (assoc exp :calls (map process-call (:calls exp))))]
    (assoc config :experiments
           (map process-exp (:experiments config)))))

(defn- no-duplicate-names?
  "Do not allow duplicate names in experiments."
  [config]
  (letfn [(exp-no-duplicate? [exp]
            (every? (fn [[_ x]] (= 1 x)) (frequencies (map :name (:calls exp)))))]
    (every? exp-no-duplicate? (:experiments config))))

(defn load-config
  "Load configuration file, handling conversion of relative to absolute paths."
  [config-file]
  {:post [(no-duplicate-names? %)]}
  (let [config (-> config-file slurp yaml/parse-string)
        base-dir (fs/file (get-in config [:dir :base] "."))
        to-process #{[:dir :out] [:dir :prep]
                     [:experiments :ref] [:experiments :intervals]
                     [:experiments :align] [:experiments :calls :file]
                     [:experiments :calls :align] [:experiments :calls :annotate]
                     [:experiments :calls :dir]}]
    (letfn [(make-absolute [x]
              (if (.isAbsolute (file x))
                x
                (str (fs/file base-dir x))))
            (maybe-process [val path]
              (if (contains? to-process path)
                (cond
                 (seq? val) (map make-absolute val)
                 (string? val) (make-absolute val)
                 :else val)
                val))
            (update-tree [config path]
              (cond (map? config)
                    (reduce (fn [item [k v]]
                              (assoc item k (cond
                                             (map? v) (update-tree v (conj path k))
                                             (seq? v) (map #(update-tree % (conj path k)) v)
                                             :else (maybe-process v (conj path k)))))
                            config
                            (vec config))
                    (contains? to-process path) (maybe-process config path)
                    :else config))]
      (-> config
          (update-tree [])
          (add-dir-files #{".vcf"})
          configure-logging))))
