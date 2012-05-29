(ns bcbio.variation.config
  "Load and prepare inputs from YAML configuration files."
  (:use [clojure.java.io])
  (:require [clj-yaml.core :as yaml]
            [fs.core :as fs]))

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
          (add-dir-files #{".vcf"})))))
