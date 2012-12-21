(ns bcbio.variation.annotate.effects
  "Predict functional consequences of variant changes leveraging snpEff."
  (:import [ca.mcgill.mcb.pcingola.snpEffect.commandLine
            SnpEffCmdEff SnpEffCmdDownload])
  (:require [clojure.java.io :as io]
            [clojure.java.shell :as shell]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]))

;; ## snpEff

(defn- get-snpeff-config
  [base-dir]
  (let [data-dir (str (fs/file base-dir "snpeff" "data"))
        orig-config (-> (ClassLoader/getSystemClassLoader)
                        (.getResourceAsStream "snpEff.config"))
        config-file (str (fs/file base-dir "snpeff" "snpEff.config"))]
    (when-not (fs/exists? data-dir)
      (fs/mkdirs data-dir))
    (when (itx/needs-run? config-file)
      (with-open [rdr (io/reader orig-config)
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
        (.parseArgs (into-array ["-c" config-file genome]))
        .run)
      (doseq [x (fs/glob (str "*" genome ".zip"))]
        (fs/delete x)))
    config-file))

(defn snpeff-annotate
  "Annotate the input file with snpEff, providing predictions of variant effects. "
  [in-file genome base-dir & {:keys [out-dir]}]
  (let [config-file (download-genome genome base-dir)
        out-file (itx/add-file-part in-file "effects" out-dir)]
    (when (itx/needs-run? out-file)
      ;; snpEff prints to standard out so we need to safely redirect that to a file.
      (let [orig-out System/out]
        (try
          (itx/with-tx-file [tx-out out-file]
            (with-open [wtr (java.io.PrintStream. tx-out)]
              (System/setOut wtr)
              (doto (SnpEffCmdEff.)
                (.parseArgs (into-array ["-noStats" "-c" config-file genome in-file]))
                .run)))
          (finally
           (System/setOut orig-out)))))
    out-file))

;; ## VEP

(defn- get-vep-cmd [vep-dir]
  (let [vep-file (when vep-dir (str (fs/file (fs/expand-home vep-dir)
                                             "variant_effect_predictor.pl")))]
    (when (and vep-file (fs/exists? vep-file))
      vep-file)))

(defn run-vep
  "Run Ensembl Variant Effects Predictor on input variant file.
   Re-annotates the input file with CSQ field compatible with Gemini."
  [in-file vep-dir & {:keys [re-run?]}]
  (when-let [vep-cmd (get-vep-cmd vep-dir)]
    (let [out-file (itx/add-file-part in-file "vep")]
      (when (or (itx/needs-run? out-file) re-run?)
        (itx/with-tx-file [tx-out out-file]
          (shell/sh "perl" vep-cmd "-i" in-file "-o" tx-out "--vcf" "--cache"
                    "--terms" "so" "--sift" "b" "--polyphen" "b" "--hgnc" "--numbers"
                    "--fields" "Consequence,Codons,Amino_acids,Gene,HGNC,Feature,EXON,PolyPhen,SIFT")))
      out-file)))
