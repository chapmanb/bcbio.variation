(ns bcbio.variation.utils.quickcompare
  "Provide a faster way to compare two VCF files, making additional use of memory.
   This is useful for comparing a larger file (whole genome sequencing) to a
   smaller one (array genotyping) where the smaller positions and genotypes
   can be loaded into memory.
   This could also be split by chromosome to reduce memory requirement but
   first pass uses the entire file."
  (:require [clojure.data.csv :as csv]
            [clojure.java.io :as io]
            [clojure.tools.cli :refer [cli]]
            [me.raynes.fs :as fs]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.run.parallel :refer [rmap]]
            [bcbio.variation.variantcontext :as gvc]))

;; ## Sample names from headers

(defn- sample-name-from-vcf
  "Retrieve first sample name from VCF file."
  [vcf-file]
  (first (.getGenotypeSamples (gvc/get-vcf-header vcf-file))))

(defn- has-sample?
  [sample vcf-file]
  (contains? (set (.getGenotypeSamples (gvc/get-vcf-header vcf-file)))
             sample))

;; ## Conversion to minimum genotypes

(defn- vc->call
  "Identify genotype call for the sample in the variant context"
  [sample vc]
  (->> (.getGenotype (:vc vc) sample)
       .getAlleles
       (map #(.getDisplayString %))
       sort))

(defn- vc->min-genotype
  "Represent a single variant context as a minimal variant call with genotype"
  [sample vc]
  (-> vc
      (select-keys [:chr :start])
      (assoc :call (vc->call sample vc))))

(defn- vcf->min-genotype
  "Convert VCF file into an in-memory representation of genotypes."
  [sample vcf-file ref-file]
  (when (has-sample? sample vcf-file)
    (with-open [vcf-iter (gvc/get-vcf-iterator vcf-file ref-file)]
      (reduce (fn [coll g]
                (assoc coll [(:chr g) (:start g)] (:call g)))
              {} (map (partial vc->min-genotype sample)
                      (gvc/parse-vcf vcf-iter))))))

;; ## Comparisons

(defn- cmp-vc-to-base-call
  "Compare variant call to a set of minimum genotypes.
   TODO: hook in flexible comparisons from multisample/compare-genotypes"
  [vc sample min-genotypes]
  (when-let [cmp (get min-genotypes [(:chr vc) (:start vc)])]
    (if (= cmp (vc->call sample vc))
      :concordant
      :discordant)))

;; ## Top level access

(defn two-vcfs
  "Compare two VCFs, base VCF will be loaded into memory and used to subset cmp-vcf."
  [base-vcf cmp-vcf ref-file]
  (let [sample (sample-name-from-vcf cmp-vcf)
        base-genotypes (vcf->min-genotype sample base-vcf ref-file)]
    (with-open [vcf-iter (gvc/get-vcf-iterator cmp-vcf ref-file)]
      (reduce (fn [coll vc]
                (if-let [cmp (cmp-vc-to-base-call vc sample base-genotypes)]
                  (assoc coll cmp (inc (get coll cmp 0)))
                  coll))
              {:sample sample}
              (gvc/parse-vcf vcf-iter)))))

(defn- percent-concordant
  [cmp]
  (when-not (zero? (apply + (vals (dissoc cmp :sample))))
    (* 100.0
       (/ (:concordant cmp)
          (apply + (vals (dissoc cmp :sample)))))))

(defn vcfdir-to-base
  "Compare a directory of single-sample VCFs to a multi-sample base VCF."
  [base-vcf vcf-dir ref-file cores]
  (let [out-file (-> (fsp/add-file-part base-vcf "cmp" vcf-dir)
                     fsp/file-root
                     (str ".csv"))
        vcf-files (map str (fs/glob (fs/file vcf-dir "*.vcf")))]
    (itx/with-tx-file [tx-out-file out-file]
      (with-open [wtr (io/writer tx-out-file)]
        (csv/write-csv wtr [["sample" "concordant" "discordant" "percent"]])
        (doseq [cmp (rmap #(two-vcfs base-vcf % ref-file) vcf-files cores)]
          (csv/write-csv wtr [[(:sample cmp) (get cmp :concordant 0) (get cmp :discordant 0)
                               (percent-concordant cmp)]]))))
    out-file))

(defn -main [& args]
  (let [[options args banner]
        (cli args
             ["-c" "--cores" "Number of cores to use" :default 1
              :parse-fn #(Integer. %)])]
    (when (not= (count args) 3)
      (println banner)
      (println "Usage: quickcompare <base-vcf> <directory of VCFs> <genome reference>")
      (System/exit 1))
    (vcfdir-to-base (first args) (second args) (nth args 2) (:cores options))))
