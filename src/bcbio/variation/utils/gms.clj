(ns bcbio.variation.utils.gms
  "Build reference Genomic Mappability Score (GMS) variant file.
  Uses full GMS files to generate VCF of potentially problematic low-GMS regions:
  http://sourceforge.net/apps/mediawiki/gma-bio/index.php"
  (:import [org.broadinstitute.variant.variantcontext VariantContextBuilder Allele]
           [org.broadinstitute.variant.variantcontext.writer VariantContextWriterFactory]
           [org.broadinstitute.variant.vcf VCFHeader
            VCFInfoHeaderLine VCFHeaderLineCount VCFHeaderLineType])
  (:use [clojure.java.io]
        [ordered.map :only [ordered-map]]
        [bcbio.align.ref :only [get-seq-dict]]
        [bcbio.variation.combine :only [combine-variants]]
        [bcbio.variation.config :only [load-config]]
        [bcbio.variation.normalize :only [hg19-map]]
        [bcbio.variation.utils.background :only [make-work-dirs]])
  (:require [clojure.java.shell :as shell]
            [clojure.string :as string]
            [me.raynes.fs :as fs]
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
    (into (ordered-map)
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

(defn- low-gms-score? [config gms-data]
  (let [thresh (get config :max-gms-score 50.0)]
    (and (> (:score gms-data) 0.0)
         (< (:score gms-data) thresh))))

(defn- gms-scores-to-vc
  "Prepare variant context from set of GMS scores"
  [techs scores]
  (let [contig (get hg19-map (-> scores first :chrom))
        all-pos (filter pos? (map :pos scores))
        pos (if (= 1 (count (set all-pos)))
              (first all-pos)
              (throw (Exception. (str "Multiple positions found: " all-pos))))
        base (->> (map :base scores)
                  (filter #(not= % "*"))
                  first)]
    (when-not (or (zero? pos) (nil? base))
      (-> (VariantContextBuilder. contig contig pos pos [(Allele/create base true)])
          (.attributes (reduce (fn [coll [tech score]]
                                 (assoc coll (str "GMS_" tech) (format "%.1f" (:score score))))
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
        (with-open [writer (VariantContextWriterFactory/create (file out-file)
                                                               (get-seq-dict ref))]
          (.writeHeader writer (get-vcf-header (keys gms-files)))
          (loop [line-iters (->> (map line-seq readers)
                                 (map (fn [x] (drop-while #(.startsWith % "#") x))))]
            (when-not (or (empty? (first line-iters))
                          (some nil? (map first line-iters)))
              (let [cur-gms (map (comp parse-gms-line first) line-iters)]
                (when (some (partial low-gms-score? ftp-config) cur-gms)
                  (when-let [vc (gms-scores-to-vc (keys gms-files) cur-gms)]
                    (.add writer vc))))
              (recur (map rest line-iters))))
          (doseq [x readers]
            (.close x)))
        (doseq [x (vals gms-files)]
          (itx/remove-path x))))
    (str out-file)))

(defn prepare-gms-vcfs
  "Prepare individual chromosome VCF files with low GMS data by sequencing technology."
  [config]
  (let [ref (:ref config)
        out-dir (get-in config [:dir :out])
        gms-by-chrom (doall (map #(prepare-vcf-at-chrom % (:ftp config) ref out-dir)
                                  (get-in config [:ftp :chromosomes])))]
    (println gms-by-chrom)
    (combine-variants gms-by-chrom ref :merge-type :full :out-dir out-dir
                      :quiet-out? true))
  (shutdown-agents))

(defn prepare-gms-vcfs-from-config [config-file]
  (let [config (load-config config-file)]
    (make-work-dirs config)
    (prepare-gms-vcfs config)))
