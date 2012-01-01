;; Functionality for running idempotent, transactional processes
;; Provides an API for long running processes in computational
;; pipelines. Avoids re-running a process if it has produced the
;; output file on a previous run, and leaving partially finished
;; files in the case of premature termination.

(ns bcbio.run.itx
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
