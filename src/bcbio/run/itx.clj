
(ns bcbio.run.itx
   "Functionality for running idempotent, transactional processes.
   Provides an API for long running processes in computational
   pipelines. Avoids re-running a process if it has produced the
   output file on a previous run, and leaving partially finished
   files in the case of premature termination."
  (:import (java.io File))
  (:use [clojure.java.io])
  (:require [fs.core :as fs]))

;; ## Idempotent processing
;; avoid re-running when output files exist

(defn needs-run?
  "Check if an output files need a run: any do not exist or empty file"
  [& fnames]
  (letfn [(file-non-empty? [f]
            (and (fs/exists? f)
                 (> (fs/size f) 0)))]
    (not-every? true?
                (map file-non-empty? (flatten fnames)))))

(defn subs-kw-files
  "Substitute any keywords in the arguments from file information map."
  [args file-info]
  (letfn [(maybe-sub-kw [x]
            (if (and (keyword? x)
                     (contains? file-info x))
              (get file-info x)
              x))]
    (map maybe-sub-kw args)))

;; ## Transactions
;; Handle output files in a separate transaction directory to avoid
;; partially finished output files if long-running processes fail.

(defn temp-dir-w-prefix [root prefix]
  (let [dir (File/createTempFile prefix "" (file root))]
    (fs/delete dir)
    (fs/mkdir dir)
    dir))

(defn safe-tx-files
  "Update file-info with need-tx files in a safe transaction directory."
  [file-info need-tx]
  (let [tx-files (map #(get file-info %) need-tx)
        tx-dir (temp-dir-w-prefix (fs/parent (first tx-files)) "txtmp")]
    (reduce (fn [m [k v]]
              (assoc m k v))
            file-info
            (zipmap need-tx
                    (map #(str (fs/file tx-dir (fs/base-name %))) tx-files)))))

(defmacro with-tx-files
  "Perform action with files, keeping need-tx files in a transaction."
  [[tx-file-info file-info need-tx] & body]
  (if (= (count need-tx) 0)
    `(do ~@body)
    `(let [~tx-file-info (safe-tx-files ~file-info ~need-tx)]
       (try
         ~@body
         (doseq [tx-key# ~need-tx]
           (fs/rename (get ~tx-file-info tx-key#) (get ~file-info tx-key#)))
         (finally
          (fs/delete-dir (fs/parent (get ~tx-file-info (first ~need-tx)))))))))

;; ## Naming
;; Generate new file names from existing ones

(defn file-root
  "Retrieve file name without extension: /path/to/fname.txt -> /path/to/fname"
  [fname]
  (let [i (.lastIndexOf fname ".")]
    (if (pos? i)
      (subs fname 0 i)
      fname)))

(defn add-file-part
  "Add file extender: base.txt -> base-part.txt"
  ([fname part]
     (add-file-part fname part nil))
  ([fname part out-dir]
     (let [out-fname (format "%s-%s%s" (file-root fname) part (fs/extension fname))]
       (if-not (nil? out-dir)
         (str (fs/file out-dir (fs/base-name out-fname)))
         out-fname))))
