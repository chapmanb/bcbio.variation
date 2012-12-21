(ns bcbio.variation.annotate.effects
  "Predict functional consequences of variant changes leveraging snpEff."
  (:import [ca.mcgill.mcb.pcingola.snpEffect.commandLine
            SnpEffCmdEff SnpEffCmdDownload])
  (:require [clojure.java.io :as io]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]))

(defn- get-snpeff-config
  [base-dir]
  (let [data-dir (str (fs/file base-dir "snpeff" "data"))
        orig-config-file (-> (ClassLoader/getSystemClassLoader)
                             (.getResource "snpEff.config")
                             (.getFile)
                             str)
        config-file (str (fs/file base-dir "snpeff" "snpEff.config"))]
    (when-not (fs/exists? data-dir)
      (fs/mkdirs data-dir))
    (when (or (itx/needs-run? config-file)
              (> (fs/mod-time orig-config-file) (fs/mod-time config-file)))
      (with-open [rdr (io/reader orig-config-file)
                  wtr (io/writer config-file)]
        (doseq [line (line-seq rdr)]
          (.write wtr (str
                       (if (.startsWith line "data_dir")
                         (str "data_dir = " data-dir)
                         line)
                       "\n")))))
    {:data-dir data-dir
     :config-file config-file}))

(defn download-genome
  "Check for a snpEff genome index and download if not present."
  [genome base-dir]
  (let [{:keys [data-dir config-file]} (get-snpeff-config base-dir)
        genome-dir (str (fs/file data-dir genome))]
    (when-not (fs/exists? genome-dir)
      (fs/mkdirs genome-dir)
      (doto (SnpEffCmdDownload.)
        (.parseArgs (into-array ["-noStats" "-c" config-file genome]))
        .run)
      (doseq [x (fs/glob (str "*" genome ".zip"))]
        (fs/delete x)))
    config-file))

(defn snpeff-annotation
  "Annotate the input file with snpEff, providing predictions of variant effects."
  [in-file genome base-dir]
  (let [config-file (download-genome genome base-dir)
        out-file (itx/add-file-part in-file "effects")]
    (when (itx/needs-run? out-file)
      (with-open [wtr (io/writer out-file)]
        (binding [*out* wtr]
          (doto (SnpEffCmdEff.)
            (.parseArgs (into-array ["-c" config-file genome in-file]))
            .run))))
    out-file))