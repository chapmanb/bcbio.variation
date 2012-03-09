(ns bcbio.variation.normalize
  "Prepare a VCF file for comparison by normalizing chromosome names,
  sort order, sample name, and genotype representation.
  This handles the work of making slightly different representations
  match, enabling VCF comparisons.
  Currently implemented for human only, with hooks to generalize for other
  organisms."
  (:import [org.broadinstitute.sting.utils.variantcontext VariantContextBuilder Genotype]
           [org.broadinstitute.sting.utils.codecs.vcf VCFHeader]
           [org.broad.tribble.readers AsciiLineReader])
  (:use [clojure.java.io]
        [bcbio.variation.variantcontext :only (write-vcf-w-template
                                               get-seq-dict get-vcf-source
                                               get-vcf-retriever
                                               get-vcf-line-parser)]
        [ordered.map :only (ordered-map)]
        [ordered.set :only (ordered-set)])
  (:require [clojure.string :as string]
            [bcbio.run.itx :as itx]
            [fs.core :as fs]))

;; ## Chromosome name remapping

;; Provide mapping from variant chromosome names to reference
;; keyed on the organism name. Currently only a human GRCh37 remap.

(defmulti chr-name-remap (fn [type & args] type))

(comment
(defn- get-hg19-map
  "Function to retrieve hg19 information. Requires korma and
  mysql connector."
  []
  (defdb db (mysql {:db "hg19"
                    :user "genome"
                    :host "genome-mysql.cse.ucsc.edu"}))
  (defentity ucscToEnsembl)
  (->> (select ucscToEnsembl)
       (map (juxt :ucsc :ensembl))
       (into {})))
)

;; Cached version of hg19 map to avoid having to make database connections
(def hg19-map
  {"chrM" "MT", "chrUn_gl000211" "GL000211", "chrUn_gl000222" "GL000222",
   "chrUn_gl000233" "GL000233", "chrUn_gl000244" "GL000244", "chrUn_gl000212" "GL000212",
   "chrUn_gl000223" "GL000223", "chrUn_gl000234" "GL000234", "chrUn_gl000245" "GL000245",
   "chrUn_gl000213" "GL000213", "chrUn_gl000224" "GL000224", "chrUn_gl000235" "GL000235",
   "chrUn_gl000246" "GL000246", "chr6_mcf_hap5" "HSCHR6_MHC_MCF", "chrUn_gl000214" "GL000214",
   "chrUn_gl000225" "GL000225", "chrUn_gl000236" "GL000236", "chrUn_gl000247" "GL000247",
   "chr1" "1", "chr6_cox_hap2" "HSCHR6_MHC_COX", "chrUn_gl000215" "GL000215",
   "chrUn_gl000226" "GL000226", "chrUn_gl000237" "GL000237", "chrUn_gl000248" "GL000248",
   "chr2" "2", "chrUn_gl000216" "GL000216", "chrUn_gl000227" "GL000227",
   "chrUn_gl000238" "GL000238", "chrUn_gl000249" "GL000249", "chr3" "3",
   "chrUn_gl000217" "GL000217", "chrUn_gl000228" "GL000228", "chrUn_gl000239" "GL000239",
   "chr9_gl000201_random" "GL000201", "chr4" "4", "chr11_gl000202_random" "GL000202",
   "chrUn_gl000218" "GL000218", "chrUn_gl000229" "GL000229", "chr9_gl000200_random" "GL000200",
   "chr19_gl000209_random" "GL000209", "chr5" "5", "chrUn_gl000219" "GL000219",
   "chr1_gl000192_random" "GL000192", "chr18_gl000207_random" "GL000207", "chr6" "6",
   "chr21_gl000210_random" "GL000210", "chr17_gl000206_random" "GL000206",
   "chr9_gl000199_random" "GL000199", "chr1_gl000191_random" "GL000191",
   "chr4_gl000194_random" "GL000194", "chr19_gl000208_random" "GL000208",
   "chr17_gl000205_random" "GL000205", "chr7" "7", "chr9_gl000198_random" "GL000198",
   "chr8_gl000197_random" "GL000197", "chr4_gl000193_random" "GL000193",
   "chr17_gl000204_random" "GL000204", "chr8" "8", "chrX" "X", "chr8_gl000196_random" "GL000196",
   "chr7_gl000195_random" "GL000195", "chr20" "20", "chr9" "9", "chrY" "Y",
   "chr17_gl000203_random" "GL000203", "chr10" "10", "chr21" "21", "chr6_dbb_hap3" "HSCHR6_MHC_DBB",
   "chr11" "11", "chr22" "22", "chr6_ssto_hap7" "HSCHR6_MHC_SSTO", "chr17_ctg5_hap1" "HSCHR17_1",
   "chr12" "12", "chr13" "13", "chr14" "14", "chr15" "15", "chr16" "16",
   "chr6_mann_hap4" "HSCHR6_MHC_MANN", "chr17" "17", "chr18" "18", "chr19" "19",
   "chr6_qbl_hap6" "HSCHR6_MHC_QBL", "chr6_apd_hap1" "HSCHR6_MHC_APD",
   "chrUn_gl000240" "GL000240", "chrUn_gl000230" "GL000230", "chrUn_gl000241" "GL000241",
   "chr4_ctg9_hap1" "HSCHR4_1", "chrUn_gl000220" "GL000220", "chrUn_gl000231" "GL000231",
   "chrUn_gl000242" "GL000242", "chrUn_gl000221" "GL000221", "chrUn_gl000232" "GL000232",
   "chrUn_gl000243" "GL000243"})

(defmethod chr-name-remap :GRCh37
  [_ ref-chrs vcf-chrs]
  (letfn [(maybe-remap-name [x]
            {:post [(contains? ref-chrs %)]}
            (if (contains? ref-chrs x)
              x
              (get hg19-map x)))]
    (zipmap vcf-chrs
            (map maybe-remap-name vcf-chrs))))

;; ## Normalize variant contexts

(defn- fix-vc
  "Build a new variant context with updated sample name."
  [sample orig]
  (letfn [(update-genotype-sample [vc]
            {:pre [(= 1 (count (.getGenotypes vc)))]}
            (let [g (first (.getGenotypes vc))]
              [(Genotype/modifyName g sample)]))]
    (-> orig
        (assoc :vc
          (-> (VariantContextBuilder. (:vc orig))
              (.genotypes (update-genotype-sample (:vc orig)))
              .make)))))

(defn- no-call-genotype?
  "Check if a variant has a non-informative no-call genotype."
  [vc]
  {:pre [(= 1 (count (:genotypes vc)))]}
  (contains? #{"NO_CALL" "MIXED" "HOM_REF"}
             (-> vc :genotypes first :type)))

(defn- sort-by-position
  "Sort stream of line inputs by position.
  Requires loading the entire file into memory during the sort-by phase
  so will not work on massive files. Should be feasible with files
  split by chromosome."
  [line-seq]
  (letfn [(add-position [line]
            [(vec (take 2 (string/split line #"\t"))) line])]
    (->> line-seq
         (map add-position)
         (sort-by first)
         (map second))))

(defn- ordered-vc-iter
  "Provide VariantContexts ordered by chromosome and normalized."
  [rdr vcf-decoder sample config]
  (->> rdr
       line-seq
       (#(if (:sort-pos config) (sort-by-position %) %))
       (map vcf-decoder)
       (remove no-call-genotype?)
       (map (partial fix-vc sample))
       (map :vc)))

(defn- fix-vcf-line
  "Provide fixes to VCF input lines that do not require VariantContext parsing.
  Fixes:
    - INFO lines with empty attributes (starting with ';'), found in
      Complete Genomics VCF files
    - Chromosome renaming."
  [line ref-info config]
  (letfn [(empty-attribute-info [info]
            (if (.startsWith info ";")
              (subs info 1)
              info))
          (fix-info [xs]
            (assoc xs 7 (-> (nth xs 7)
                            empty-attribute-info)))
          (fix-chrom [new xs]
            (assoc xs 0 new))]
    (let [parts (string/split line #"\t")
          cur-chrom (first (vals
                            (chr-name-remap (:org config) ref-info [(first parts)])))]
      {:chrom cur-chrom
       :line (->> parts
                  (fix-chrom cur-chrom)
                  fix-info
                  (string/join "\t"))})))

(defn- vcf-by-chrom
  "Split input VCF into separate files by chromosome, returning a map of file names."
  [vcf-file ref-file tmp-dir config]
  (letfn [(ref-chr-files [ref-file]
            (into (ordered-map)
                  (map (fn [x] [(.getSequenceName x)
                                (str (fs/file tmp-dir (str "prep" (.getSequenceName x) ".vcf")))])
                       (-> ref-file get-seq-dict .getSequences))))
          (write-by-chrom [ref-wrtrs line]
            (let [line-info (fix-vcf-line line ref-wrtrs config)]
              (.write (get ref-wrtrs (:chrom line-info))
                      (str (:line line-info) "\n"))))]
    (let [ref-chrs (ref-chr-files ref-file)
          ref-wrtrs (zipmap (keys ref-chrs) (map writer (vals ref-chrs)))]
      (with-open [rdr (reader vcf-file)]
        (itx/with-open-map ref-wrtrs
          (->> rdr
               line-seq
               (drop-while #(.startsWith % "#"))
               (map (partial write-by-chrom ref-wrtrs))
               doall))
        ref-chrs))))

;; ## Top level functionality to manage inputs and writing.

(defn- update-header
  "Update header information, removing contig and adding sample names."
  [sample]
  (fn [header]
    {:pre [(= 1 (count (.getGenotypeSamples header)))]}
    (VCFHeader. (apply ordered-set (remove #(= "contig" (.getKey %)) (.getMetaData header)))
                (ordered-set sample))))

(defn- write-prepped-vcf
  "Write VCF file with correctly ordered and cleaned variants."
  [vcf-file out-info ref-file sample config]
  (itx/with-temp-dir [tmp-dir (fs/parent (:out out-info))]
    (let [reader-by-chr (into (ordered-map) (map (fn [[k v]] [k (reader v)])
                                                 (vcf-by-chrom vcf-file ref-file tmp-dir config)))]
      (itx/with-open-map reader-by-chr
        (with-open [vcf-reader (AsciiLineReader. (input-stream vcf-file))]
          (let [vcf-decoder (get-vcf-line-parser vcf-reader)]
            (write-vcf-w-template vcf-file out-info
                                  (flatten
                                   (for [rdr (vals reader-by-chr)]
                                     (ordered-vc-iter rdr vcf-decoder sample config)))
                                  ref-file
                                  :header-update-fn (update-header sample))))))))

(defn prep-vcf
  "Prepare VCF for comparison by normalizing high level attributes
  Assumes by position sorting of variants in the input VCF. Chromosomes do
  not require a specific order, but positions internal to a chromosome do.
  Currently configured for human preparation."
  [in-vcf-file ref-file sample & {:keys [out-dir out-fname]}]
  (let [config {:org :GRCh37 :sort-pos true}
        base-name (if (nil? out-fname) (itx/remove-zip-ext in-vcf-file) out-fname)
        out-file (itx/add-file-part base-name "prep" out-dir)]
    (if (itx/needs-run? out-file)
      (write-prepped-vcf in-vcf-file {:out out-file} ref-file sample config))
    out-file))
