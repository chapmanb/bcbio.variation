(ns bcbio.variation.utils.background
  "Prepare annotated VCF files to use as background for variant calling and recalibration.
  Batch variant calling and recalibration with GATK improves resulting calls. This
  provides a ready to use set of calls to batch with a single sample using 1000 genomes data."
  (:use [clojure.java.io]
        [bcbio.variation.compare :only [load-config]]
        [bcbio.variation.combine :only [select-by-sample combine-variants]])
  (:require [clojure.java.shell :as shell]
            [clj-yaml.core :as yaml]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]))

(defn- download-sample-bam
  "Download BAM file and index for a sample from 1000 genomes FTP."
  [sample ftp-config out-dir]
  (letfn [(download [url fname]
            (if-not (fs/exists? fname)
              (println "Downloading" url "to" fname)
              (shell/with-sh-dir out-dir
                (shell/sh "wget" "-O" fname url))))]
    (let [dl-url (format (:bam-url ftp-config) sample sample)
          local-file (str (fs/file out-dir (fs/base-name dl-url)))]
      (download dl-url local-file)
      (download (str dl-url ".bai") (str local-file ".bai"))
      local-file)))

(defn- merge-sample-batch
  "Select sample from input VCF batch and merge."
  [sample vcfs ref out-dir]
  (combine-variants (map #(select-by-sample sample (:file %) (:chrom %)
                                            ref :out-dir out-dir)
                                   vcfs)
                      ref :out-dir out-dir))

(defn- download-base-vcfs
  "Download original VCFs from 1000 genomes for processing."
  [chroms ftp-config out-dir]
  (letfn [(download-vcf [url fname]
            (println "Downloading" url "to" fname)
            (shell/with-sh-dir out-dir
              (shell/sh "wget" "-O" fname url)
              (shell/sh "gunzip" (str (fs/base-name fname)))))]
    (for [chrom chroms]
      (let [dl-url (format (:vcf-url ftp-config) chrom)
            local-file (str (fs/file out-dir (fs/base-name dl-url)))
            final-file (itx/remove-zip-ext local-file)]
        (when-not (fs/exists? final-file)
          (download-vcf dl-url local-file))
        {:chrom chrom :file final-file}))))

(defn make-work-dirs [config]
  (doseq [dir-name (-> config :dir keys)]
    (let [cur-dir (get-in config [:dir dir-name])]
      (when-not (fs/exists? cur-dir)
        (fs/mkdirs cur-dir)))))

(defn -main [config-file]
  (let [config (load-config config-file)]
    (make-work-dirs config)
    (doall
     (for [chroms (partition-all (get-in config [:ftp :batch-size])
                                 (get-in config [:ftp :chromosomes]))]
       (let [vcfs (download-base-vcfs chroms (:ftp config)
                                      (get-in config [:dir :prep]))]
         (for [sample (:genomes config)]
           (merge-sample-batch sample vcfs (:ref config)
                                  (get-in config [:dir :prep]))))))))
