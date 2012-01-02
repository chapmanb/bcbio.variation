;; Functionality for running idempotent, transactional processes
;; Provides an API for long running processes in computational
;; pipelines. Avoids re-running a process if it has produced the
;; output file on a previous run, and leaving partially finished
;; files in the case of premature termination.

(ns bcbio.run.itx
  (:import (java.io File))
  (:use [clojure.java.io])
  (:require [fs.core :as fs]))

;; Idempotent processing: avoid re-running when output files exist

(defn needs-run? [& fnames]
  "Check if an output files need a run: any do not exist or empty file"
  (letfn [(file-non-empty? [f]
            (and (fs/exists? f)
                 (> (fs/size f) 0)))]
    (not-every? true?
                (map file-non-empty? (flatten fnames)))))

(defn subs-kw-files [args file-info]
  "Substitute any keywords in the arguments from file information map."
  (letfn [(maybe-sub-kw [x]
            (if (and (keyword? x)
                     (contains? file-info x))
              (get file-info x)
              x))]
    (map maybe-sub-kw args)))

;; Handle output files in a separate transaction directory to avoid
;; partially finished output files if long-running processes fail.

(defn temp-dir-w-prefix [root prefix]
  (let [dir (File/createTempFile prefix "" (file root))]
    (fs/delete dir)
    (fs/mkdir dir)
    dir))

(defn safe-tx-files [file-info need-tx]
  "Update file-info with need-tx files in a safe transaction directory."
  (let [tx-files (map #(get file-info %) need-tx)
        tx-dir (temp-dir-w-prefix (fs/parent (first tx-files)) "txtmp")]
    (reduce (fn [m [k v]]
              (assoc m k v))
            file-info
            (zipmap need-tx
                    (map #(str (fs/file tx-dir (fs/base-name %))) tx-files)))))

(defmacro with-tx-files [[tx-file-info file-info need-tx] & body]
  "Perform action with files, keeping need-tx files in a transaction."
  (if (= (count need-tx) 0)
    `(do ~@body)
    `(let [~tx-file-info (safe-tx-files ~file-info ~need-tx)]
       (try
         ~@body
         (doseq [tx-key# ~need-tx]
           (fs/rename (get ~tx-file-info tx-key#) (get ~file-info tx-key#)))
         (finally
          (fs/delete-dir (fs/parent (get ~tx-file-info (first ~need-tx)))))))))

;; Utility functions for generating new file names from existing ones

(defn file-root [fname]
  "Retrieve file name without extension: /path/to/fname.txt -> /path/to/fname"
  (let [i (.lastIndexOf fname ".")]
    (if (pos? i)
      (subs fname 0 i)
      fname)))

(defn add-file-part [fname part]
  "Add file extender: base.txt -> base-part.txt"
  (format "%s-%s%s" (file-root fname) part (fs/extension fname)))
