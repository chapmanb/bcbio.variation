(ns bcbio.variation.utils.gms
  "Build reference Genomic Mappability Score (GMS) variant file.
  Uses full GMS files to generate VCF of potentially problematic low-GMS regions:
  http://sourceforge.net/apps/mediawiki/gma-bio/index.php"
  (:import [org.broadinstitute.sting.utils.variantcontext VariantContextBuilder Allele]
           [org.broadinstitute.sting.utils.codecs.vcf StandardVCFWriter VCFHeader
            VCFInfoHeaderLine VCFHeaderLineCount VCFHeaderLineType])
  (:use [clojure.java.io]
        [bcbio.align.ref :only [get-seq-dict]]
        [bcbio.variation.config :only [load-config]]
        [bcbio.variation.normalize :only [hg19-map]]
        [bcbio.variation.utils.background :only [make-work-dirs]])
  (:require [clojure.java.shell :as shell]
            [clojure.string :as string]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]))

(defn- download-chrom-gms-data
  "Download GMS data for all technologies at a chromosome."
  [chrom ftp-config out-dir]
  (letfn [(download-gms [chrom tech]
            (let [dl-url (format (:gms-url ftp-config) (:genome-build ftp-config)
                                 tech chrom)
                  final-file (itx/add-file-part (itx/remove-zip-ext (fs/base-name dl-url))
                                                tech out-dir)
                  dl-file (str final-file ".gz")]
              (when (itx/needs-run? final-file)
                (shell/with-sh-dir out-dir
                  (println (format "Downloading %s to %s" dl-url dl-file))
                  (shell/sh "wget" "-O" dl-file dl-url)
                  (shell/sh "gunzip" dl-file)))
              final-file))]
    (into {}
          (map (juxt identity (partial download-gms chrom))
               (:technologies ftp-config)))))

(defn- parse-gms-line
  "Retrieve chromosome, position and GMS score for line in a GMS file"
  [line]
  (let [[chrom pos base _ _ score] (string/split line #"\t")]
    {:chrom chrom
     :pos (Integer/parseInt pos)
     :base base
     :score (Float/parseFloat score)}))

(defn- low-gms-score? [gms-data]
  (and (> (:score gms-data) 0.0)
       (< (:score gms-data) 50.0)))

(defn- gms-scores-to-vc
  "Prepare variant context from set of GMS scores"
  [techs scores]
  (let [contig (get hg19-map (-> scores first :chrom))
        pos (->> (map :pos scores)
                 (filter pos?)
                 first)
        base (->> (map :base scores)
                  (filter #(not= % "*"))
                  first)]
    (when-not (or (zero? pos) (nil? base))
      (-> (VariantContextBuilder. contig contig pos pos [(Allele/create base true)])
          (.attributes (reduce (fn [coll [tech score]]
                                 (assoc coll (str "GMS_" tech) (:score score)))
                               {} (map vector techs scores)))
          (.make)))))

(defn- get-vcf-header [techs]
  (VCFHeader. (set
               (map #(VCFInfoHeaderLine. (format "GMS_%s" %) 1
                                         VCFHeaderLineType/Float
                                         (format "Genome Mappability Score: %s" %))
                    techs))))

(defn- prepare-vcf-at-chrom
  "Prepare an output VCF of low GMS values at the provided chromosome."
  [chrom ftp-config ref out-dir]
  (let [out-file (file out-dir (format "lowgms-scores-%s.vcf" chrom))]
    (when (itx/needs-run? out-file)
      (let [gms-files (download-chrom-gms-data chrom ftp-config out-dir)
            readers (map reader (vals gms-files))]
        (with-open [writer (StandardVCFWriter. (file out-file) (get-seq-dict ref))]
          (.writeHeader writer (get-vcf-header (keys gms-files)))
          (loop [line-iters (->> (map line-seq readers)
                                 (map (fn [x] (drop-while #(.startsWith % "#") x))))]
            (when-not (empty? (first line-iters))
              (let [cur-gms (map (comp parse-gms-line first) line-iters)]
                (when (some low-gms-score? cur-gms)
                  (when-let [vc (gms-scores-to-vc (keys gms-files) cur-gms)]
                    (.add writer vc))))
              (recur (map rest line-iters))))
          (doseq [x readers]
            (.close x)))
        (doseq [x (vals gms-files)]
          (itx/remove-path x))))
    out-file))

(defn prepare-gms-vcfs
  "Prepare individual chromosome VCF files with low GMS data by sequencing technology."
  [config]
  (let [gms-by-chrom (doall (pmap #(prepare-vcf-at-chrom % (:ftp config) (:ref config)
                                                         (get-in config [:dir :out]))
                                  (get-in config [:ftp :chromosomes])))]
    (println gms-by-chrom))
  (shutdown-agents))

(defn -main [config-file]
  (let [config (load-config config-file)]
    (make-work-dirs config)
    (prepare-gms-vcfs config)))
