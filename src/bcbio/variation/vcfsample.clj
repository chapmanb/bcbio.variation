(ns bcbio.variation.vcfsample
  "Sort VCF sample columns to have a consistent order between multiple inputs.
   Variant callers order called outputs differently and this ensures they are
   consistent to feed into ensemble calling."
  (:import [htsjdk.samtools.util BlockCompressedInputStream BlockCompressedOutputStream])
  (:require [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.variation.variantcontext :as vc]
            [clojure.core.strint :refer [<<]]
            [clojure.java.io :as io]
            [clojure.string :as string]
            [me.raynes.fs :as fs]))

;; ## bgzipped and indexed files

(defn- tabix-index-vcf
  "Tabix index input VCF inside a transactional directory."
  [bgzip-file]
  (let [tabix-file (str bgzip-file ".tbi")]
    (when (or (itx/needs-run? tabix-file) (not (itx/up-to-date? tabix-file bgzip-file)))
      (itx/with-tx-file [tx-tabix-file tabix-file]
        (let [tx-bgzip-file (fsp/file-root tx-tabix-file)
              full-bgzip-file (str (fs/file bgzip-file))
              tmp-dir (str (fs/parent tx-bgzip-file))]
          (itx/check-run (<< "ln -s ~{full-bgzip-file} ~{tx-bgzip-file}"))
          (itx/check-run (<< "bcftools tabix -p vcf ~{tx-bgzip-file}")))))
    tabix-file))

(defn bgzip-index-vcf
  "Prepare a VCF file for positional query with bgzip and tabix indexing."
  [vcf-file & {:keys [remove-orig? remove-nopass? dir]}]
  (let [out-orig (str (fsp/file-root vcf-file)
                      (if remove-nopass? "-passonly" "")
                      ".vcf.gz")
        out-file (if dir (str (io/file dir (fs/base-name out-orig))) out-orig)]
    (if remove-nopass?
      (itx/run-cmd out-file "bcftools view -f 'PASS,.' ~{vcf-file} | bgzip -c > ~{out-file}")
      (itx/run-cmd out-file "bgzip -c ~{vcf-file} > ~{out-file}"))
    (when (and (not (.endsWith vcf-file ".gz")) remove-orig?)
      (fsp/remove-path vcf-file))
    (tabix-index-vcf out-file)
    out-file))

;; ## Unique file names

(defn unique-work-file
  "Create a work file with unique name in case of shared base names."
  [orig-file ext all-files work-dir]
  (let [cmp-files (remove #(= % orig-file) all-files)
        parts (reverse (string/split orig-file #"/"))
        unique-file (loop [i 1]
                      (let [cur-name (string/join "-" (reverse (take i parts)))]
                        (if (not-any? #(.endsWith % cur-name) cmp-files)
                          cur-name
                          (recur (inc i)))))]
    (fsp/add-file-part unique-file ext work-dir)))

;; ## Consistent order sorting

(defn- sort-sample-line
  "Sort samples in a VCF line using reordered indexes from calculate-reorder."
  [line reorder]
  (let [[keep samples] (split-at 9 (string/split line #"\t"))]
    (string/join "\t"
                 (concat keep
                         (->> samples
                              (map-indexed vector)
                              (sort-by (fn [[i x]] (get reorder i)))
                              (map second))))))

(defn- calculate-reorder
  "Create a dictionary of sample indexes in the original VCF to those desired in the output."
  [orig-order want-order]
  (let [want-indexes (reduce (fn [coll [i x]]
                               (assoc coll x i))
                             {} (map-indexed vector want-order))]
    (reduce (fn [coll [i x]]
              (assoc coll i (get want-indexes x)))
            {} (map-indexed vector orig-order))))

(defn- sort-samples
  "Sort samples in a VCF file, moving from orig-order to want-order."
  [vcf-file orig-order want-order all-vcfs work-dir]
  (let [out-file (unique-work-file vcf-file "ssort" all-vcfs work-dir)
        sample-reorder (calculate-reorder orig-order want-order)]
    (with-open [rdr (io/reader (BlockCompressedInputStream. (io/file vcf-file)))
                wtr (io/writer (BlockCompressedOutputStream. (io/file out-file)))]
      (doseq [line (line-seq rdr)]
        (.write wtr (str (if (.startsWith line "##")
                           line
                           (sort-sample-line line sample-reorder))
                         "\n"))))
    (bgzip-index-vcf out-file)
    out-file))

(defn- maybe-sort-names
  "Extract sample names for the current file and do sorting if needed."
  [vcf-file sorder all-vcfs work-dir]
  (let [cur-sorder (vc/get-samples vcf-file)]
    (if (not= cur-sorder sorder)
      (sort-samples vcf-file cur-sorder sorder all-vcfs work-dir)
      vcf-file)))

(defn consistent-order
  "Ensure the set of VCF files have a consistent sample order relative to the first."
  [vcf-files work-dir]
  (fsp/safe-mkdir work-dir)
  (let [sorder (vc/get-samples (first vcf-files))]
    (cons (first vcf-files)
          (map #(maybe-sort-names % sorder vcf-files work-dir) (rest vcf-files)))))
