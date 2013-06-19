(ns bcbio.variation.utils.background
  "Prepare annotated VCF files to use as background for variant calling and recalibration.
  Batch variant calling and recalibration with GATK improves resulting calls. This
  provides a ready to use set of calls to batch with a single sample using 1000 genomes data."
  (:use [clojure.java.io]
        [bcbio.variation.config :only [load-config]]
        [bcbio.variation.combine :only [combine-variants]]
        [bcbio.variation.filter.intervals :only [select-by-sample]]
        [bcbio.variation.annotation :only [add-gatk-annotations]])
  (:require [clojure.java.shell :as shell]
            [clojure.string :as string]
            [me.raynes.fs :as fs]
            [aws.sdk.s3 :as s3]
            [bcbio.run.itx :as itx]))

;; ## Combine and annotate VCFs

(defn- download-sample-bam
  "Download BAM file and index for a sample from 1000 genomes FTP."
  [sample ftp-config out-dir]
  (letfn [(download [url fname]
            (when-not (fs/exists? fname)
              (println "Downloading" url "to" fname)
              (shell/with-sh-dir out-dir
                (shell/sh "wget" "-O" fname url))))]
    (let [dl-url (format (:bam-url ftp-config) sample sample)
          local-file (str (fs/file out-dir (string/replace (fs/base-name dl-url) ".*" "")))]
      (download dl-url local-file)
      (download (str dl-url ".bai") (str local-file ".bai"))
      local-file)))

(defn- gzip-needs-run?
  "Check if a file exists, also checking for gzipped versions."
  [x]
  (every? itx/needs-run? [x (str x ".gz")]))

(defn- annotate-sample
  "Annotate genome sample VCFs with GATK metrics."
  [sample-info ref ftp-config prep-dir out-dir]
  (let [final-file (str (fs/file out-dir (format "%s-annotated.vcf" (:sample sample-info))))]
    (when (gzip-needs-run? final-file)
      (let [sample-bam (download-sample-bam (:sample sample-info) ftp-config prep-dir)
            ann-vcf (add-gatk-annotations (:file sample-info) sample-bam ref)]
        (fs/rename ann-vcf final-file)
        (itx/remove-path sample-bam)))
    final-file))

(defn- combine-samples
  "Combine sample VCFs split by chromosome."
  [sample-info ref out-dir]
  (letfn [(combine-sample [[name xs]]
            {:sample name
             :file (combine-variants (map :file xs) ref :merge-type :full
                                     :out-dir out-dir)})]
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

;; ## Create combined background file

(defn prep-combined-background
  "Prepare combined VCF file with background information from multiple inputs."
  [vcfs config]
  (letfn [(maybe-bgzip-vcf [x]
            (first (filter fs/exists? [x (str x ".gz")])))]
    (let [out-dir (get-in config [:dir :out])
          out-file (str (fs/file out-dir (get-in config [:upload :combined-vcf])))]
      (when (gzip-needs-run? out-file)
        (-> (combine-variants (map maybe-bgzip-vcf vcfs) (:ref config)
                              :merge-type :full :out-dir out-dir)
            (fs/rename out-file)))
      out-file)))

;; ## Tabix prep and upload

(defn- tabix-prep-vcf
  "Prep VCF for tabix access by bgzipping and indexing."
  [vcf]
  (let [out-dir (str (fs/parent vcf))
        vcf-gz (str vcf ".gz")
        tbi-gz (str vcf-gz ".tbi")]
    (shell/with-sh-dir out-dir
      (when (itx/needs-run? vcf-gz)
        (shell/sh "bgzip" (str (fs/base-name vcf))))
      (when (itx/needs-run? tbi-gz)
        (shell/sh "tabix" "-p" "vcf" (str (fs/base-name vcf-gz)))))
    [vcf-gz tbi-gz]))

(defmulti upload-result-vcf
  "Upload prepared sample VCF bgzipped and tabix indexed."
  (fn [_ config] (keyword (get-in config [:upload :target]))))

(defmethod upload-result-vcf :s3
  [vcf config]
  (let [cred {:access-key (System/getenv "AWS_ACCESS_KEY_ID")
              :secret-key (System/getenv "AWS_SECRET_ACCESS_KEY")}
        bucket (get-in config [:upload :bucket])]
    (when-not (s3/bucket-exists? cred bucket)
      (s3/create-bucket cred bucket)
      (s3/update-bucket-acl cred bucket (s3/grant :all-users :read)))
    (doseq [fname (tabix-prep-vcf vcf)]
      (let [s3-key (format "%s/%s" (get-in config [:upload :folder])
                           (str (fs/base-name fname)))]
        (when-not (s3/object-exists? cred bucket s3-key)
          (s3/put-object cred bucket s3-key (file fname))
          (s3/update-object-acl cred bucket s3-key
                                (s3/grant :all-users :read)))
        (println s3-key)))))

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
                           (sort-by :sample combo-samples))]
      (doseq [ready-vcf (cons (prep-combined-background ann-samples config)
                              ann-samples)]
        (upload-result-vcf ready-vcf config)))))
