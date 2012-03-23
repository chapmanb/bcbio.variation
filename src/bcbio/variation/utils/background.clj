(ns bcbio.variation.utils.background
  "Prepare annotated VCF files to use as background for variant calling and recalibration.
  Batch variant calling and recalibration with GATK improves resulting calls. This
  provides a ready to use set of calls to batch with a single sample using 1000 genomes data."
  (:use [clojure.java.io]
        [bcbio.variation.compare :only [load-config]])
  (:require [clojure.java.shell :as shell]
            [clj-yaml.core :as yaml]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]))

(defn- download-base-vcfs
  "Download original VCFs from 1000 genomes for processing."
  [ftp-config out-dir]
  (letfn [(download-vcf [url fname]
            (println "Downloading" url "to" fname)
            (shell/with-sh-dir out-dir
              (shell/sh "wget" "-O" fname url)
              (shell/sh "gunzip" (str (fs/base-name fname)))))]
    (doall
     (for [chrom (:chromosomes ftp-config)]
       (let [dl-url (format (:vcf_url ftp-config) chrom)
             local-file (str (fs/file out-dir (fs/base-name dl-url)))
             final-file (itx/remove-zip-ext local-file)]
         (when-not (fs/exists? final-file)
           (download-vcf dl-url local-file))
         final-file)))))

(defn make-work-dirs [config]
  (doseq [dir-name (-> config :dir keys)]
    (let [cur-dir (get-in config [:dir dir-name])]
      (when-not (fs/exists? cur-dir)
        (fs/mkdirs cur-dir)))))

(defn -main [config-file]
  (let [config (load-config config-file)]
    (make-work-dirs config)
    (let [vcfs (download-base-vcfs (:ftp config) (get-in config [:dir :prep]))]
      (println vcfs))))
