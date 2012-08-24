(ns bcbio.variation.index.gemini
  "Index and retrieve variant associated population genetic and disease data.
   Built on the Gemini framework: https://github.com/arq5x/gemini"
  (:require [clojure.java.jdbc :as sql]
            [clojure.java.shell :as shell]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]))

(defn gemini-installed? []
  (let [info (try (shell/sh "gemini" "-h")
                  (catch java.io.IOException _
                    {:exit -1}))]
    (zero? (:exit info))))

(defn index-variant-file
  "Pre-index a variant file with gemini"
  [in-file ref-file & {:keys [re-index?]}]
  (when (gemini-installed?)
    (let [index-file (str (itx/file-root in-file) "-gemini.db")]
      (when (or (itx/needs-run? index-file) re-index?)
        (itx/with-tx-file [tx-index index-file]
          (shell/sh "gemini" "load" "-v" in-file tx-index)))
      index-file)))