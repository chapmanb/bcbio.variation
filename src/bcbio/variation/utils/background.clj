(ns bcbio.variation.utils.background
  "Prepare annotated VCF files to use as background for variant calling and recalibration.
  Batch variant calling and recalibration with GATK improves resulting calls. This
  provides a ready to use set of calls to batch with a single sample using 1000 genomes data."
  (:use [clojure.java.io]
        [bcbio.variation.compare :only [load-config]]
        [bcbio.variation.combine :only [select-by-sample combine-variants]]
        [bcbio.variation.annotation :only [add-variant-annotations]])
  (:require [clojure.java.shell :as shell]
            [clj-yaml.core :as yaml]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]))

;; ## Combine and annotate VCFs

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

(defn- annotate-sample
  "Annotate genome sample VCFs with GATK metrics."
  [sample-info ref ftp-config prep-dir out-dir]
  (let [final-file (str (fs/file out-dir (format "%s-annotated.vcf" (:sample sample-info))))]
    (when (itx/needs-run? final-file)
      (let [sample-bam (download-sample-bam (:sample sample-info) ftp-config prep-dir)
            ann-vcf (add-variant-annotations (:file sample-info) sample-bam ref)]
        (fs/rename ann-vcf final-file)
        (itx/remove-path sample-bam)))
    final-file))

(defn- combine-samples
  "Combine sample VCFs split by chromosome."
  [sample-info ref out-dir]
  (letfn [(combine-sample [[name xs]]
            {:sample name
             :file (combine-variants (map :file xs) ref :out-dir out-dir)})]
    (map combine-sample
         (group-by :sample (flatten sample-info)))))

;; ## Subset VCF by sample

(defn- download-chrom-vcf
  "Download chromosome VCF from 1000 genomes for processing."
  [chrom ftp-config out-dir]
  (letfn [(download-vcf [url fname]
            (println "Downloading" url "to" fname)
            (shell/with-sh-dir out-dir
              (shell/sh "wget" "-O" fname url)
              (shell/sh "gunzip" (str (fs/base-name fname)))))]
    (let [dl-url (format (:vcf-url ftp-config) chrom)
          local-file (str (fs/file out-dir (fs/base-name dl-url)))
          final-file (itx/remove-zip-ext local-file)]
      (when-not (fs/exists? final-file)
        (download-vcf dl-url local-file))
      final-file)))

(defn- select-samples-at-chrom
  "Select samples from input 1000 genomes chromosome VCF."
  [chrom samples ref ftp-config out-dir]
  (let [sample-info (map (fn [x] {:sample x
                                  :file (str (fs/file out-dir (format "%s-%s.vcf" x chrom)))})
                         samples)]
    (when (apply itx/needs-run? (map :file sample-info))
      (let [chrom-vcf (download-chrom-vcf chrom ftp-config out-dir)]
        (doseq [sample samples]
          (select-by-sample sample chrom-vcf chrom ref :out-dir out-dir
                            :remove-refcalls true))
        (itx/remove-path chrom-vcf)))
    sample-info))

(defn make-work-dirs [config]
  (doseq [dir-name (-> config :dir keys)]
    (let [cur-dir (get-in config [:dir dir-name])]
      (when-not (fs/exists? cur-dir)
        (fs/mkdirs cur-dir)))))

(defn -main [config-file]
  (let [config (load-config config-file)]
    (make-work-dirs config)
    (let [prep-dir (get-in config [:dir :prep])
          samples (map #(select-samples-at-chrom % (:genomes config) (:ref config)
                                                 (:ftp config) prep-dir)
                                        (get-in config [:ftp :chromosomes]))
          combo-samples (combine-samples samples (:ref config) prep-dir)
          ann-samples (map #(annotate-sample % (:ref config) (:ftp config)
                                             prep-dir (get-in config [:dir :out]))
                           combo-samples)]
      (println ann-samples))))
