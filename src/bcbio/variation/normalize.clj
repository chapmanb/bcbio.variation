(ns bcbio.variation.normalize
  "Prepare a VCF file for comparison by normalizing chromosome names,
  sort order, sample name, and genotype representation.
  This handles the work of making slightly different representations
  match, enabling VCF comparisons.
  Currently implemented for human only, with hooks to generalize for other
  organisms."
  (:import [org.broadinstitute.variant.variantcontext VariantContextBuilder GenotypeBuilder]
           [org.broadinstitute.variant.vcf VCFHeader]
           [org.broad.tribble.readers AsciiLineReader]
           [org.apache.commons.lang CharUtils])
  (:use [clojure.java.io]
        [bcbio.variation.variantcontext :only [write-vcf-w-template
                                               get-vcf-iterator parse-vcf
                                               get-vcf-line-parser
                                               from-genotype]]
        [bcbio.align.ref :only [get-seq-dict get-seq-name-map extract-sequence]]
        [ordered.map :only (ordered-map)]
        [ordered.set :only (ordered-set)])
  (:require [clojure.string :as string]
            [me.raynes.fs :as fs]
            [lonocloud.synthread :as ->]
            [bcbio.run.itx :as itx]))

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
  {"chrM" "MT" "chrMT" "MT" "chrUn_gl000211" "GL000211", "chrUn_gl000222" "GL000222",
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

(defn- fix-non-version-names
  "Convert any non-versioned names into the representative version in ref-dict."
  [base-map ref-dict]
  (letfn [(find-best-match [x check]
            (first (filter #(.startsWith % x) check)))]
    (reduce (fn [coll [k v]]
              (assoc coll k
                     (if (contains? ref-dict v)
                       v
                       (find-best-match v (keys ref-dict)))))
            {} base-map)))

(defn- add-alt-keys
  "Add alternative key variations:
    - underscore to dash in hg19 names
    - chr added to all GRCh37 names instead of hg19 names"
  [base-map modtype]
  {:pre [(= modtype :underscore)]}
  (reduce (fn [coll [k v]]
            (-> coll
                (assoc k v)
                (assoc (string/replace k "_" "-") v)
                (assoc (str "chr" v) v)))
          {} base-map))

(defn prep-rename-map
  "Fix GRCh37/hg19 name mappings to handle common problem cases."
  [map-key ref-file]
  (let [remappers {:GRCh37 hg19-map}]
    (-> (get remappers map-key)
        (fix-non-version-names (get-seq-name-map ref-file))
        (add-alt-keys :underscore))))

(defn- chrs-from-fasta-file
  "Retrieve a list of all chromosome names from a reference FASTA file."
  [ref-file]
  (map #(.getSequenceName %) (-> ref-file get-seq-dict .getSequences)))

(defmethod chr-name-remap :GRCh37
  [map-key ref-file orig-ref-file]
  (let [rename-map (prep-rename-map map-key ref-file)
        ref-chrs (set (chrs-from-fasta-file ref-file))
        vcf-chrs (when (and orig-ref-file (not= orig-ref-file ref-file))
                   (chrs-from-fasta-file orig-ref-file))]
    (reduce (fn [coll x]
              (let [remap-x (get coll x)]
                (if (and remap-x (contains? ref-chrs remap-x))
                  (assoc coll x remap-x)
                  (assoc coll x x))))
            rename-map vcf-chrs)))

;; ## Resort and normalize variants

(defn- fix-vc
  "Build a new variant context with updated sample name and normalized alleles.
  Based on :prep-allele-count in the configuration updates haploid allele calls. This
  normalizes the representation in Mitochondrial and Y chromosomes which are
  haploid but are often represented as diploid with a single call."
  [sample config orig]
  (letfn [(update-genotype-sample [vc sample]
            (if (and (not (nil? sample))
                     (= 1 (count (.getGenotypes vc))))
              (let [g (first (.getGenotypes vc))]
                [(-> (GenotypeBuilder. g)
                     (.name sample)
                     .make)])
              (.getGenotypes vc)))
          (normalize-allele-calls [g]
            {:pre [(or (nil? (:prep-allele-count config))
                       (contains? (set [1 (:prep-allele-count config)]) (count (.getAlleles g))))]}
            (if (or (nil? (:prep-allele-count config))
                    (= (count (.getAlleles g)) (:prep-allele-count config)))
              g
              (-> (GenotypeBuilder. g)
                  (.alleles (repeat (:prep-allele-count config)
                                    (first (.getAlleles g))))
                  .make)))]
    (-> orig
        (assoc :vc
          (-> (VariantContextBuilder. (:vc orig))
              (.genotypes (map normalize-allele-calls (update-genotype-sample (:vc orig) sample)))
              .make)))))

(defn- no-call-genotype?
  "Check if a variant has a non-informative no-call genotype."
  [vc config]
  (let [to-remove (cond
                   (:prep-sv-genotype config) #{}
                   (:remove-refcalls config) #{"NO_CALL" "MIXED" "HOM_REF"}
                   :else #{"NO_CALL" "MIXED"})]
    (if-not (= 1 (:num-samples vc)) false
            (try
              (contains? to-remove (-> vc :genotypes first :type))
              (catch Exception e
                (println (:chr vc) (:start vc))
                (throw e))))))

(defn nochange-alt?
  "Check a VCF input line for identical REF and ALT calls"
  [line]
  (let [parts (string/split line #"\t")]
    (= (nth parts 3) (nth parts 4))))

(defn- sort-by-position
  "Sort stream of line inputs by position.
  Requires loading the entire file into memory during the sort-by phase
  so will not work on massive files. Should be feasible with files
  split by chromosome."
  [line-seq]
  (letfn [(add-position [line]
            (let [[chrom start] (take 2 (string/split line #"\t"))]
              [[chrom (Integer/parseInt start)] line]))]
    (->> line-seq
         (map add-position)
         (sort-by first)
         (map second))))

(defn- normalize-sv-genotype
  "Provide genotype calls for structural variants to a single ref call.
  Structural variants often don't have proper genotype references since
  individual haplotypes are not called. This makes them a single reference
  if not specified or mixed."
  [config sample orig]
  (letfn [(maybe-fix-vc [g alt-allele]
            (case (:type g)
              "MIXED" (-> (GenotypeBuilder. (:genotype g))
                          (.alleles) (remove #(.isNoCall %) (:alleles g))
                          .make)
              ("UNAVAILABLE" "NO_CALL") (-> (GenotypeBuilder. (:genotype g))
                                            (.alleles [alt-allele])
                                            .make)
              (:genotype g)))
          (ref-vc-genotype [gs alt-allele]
            (case (count gs)
              0 [(-> (GenotypeBuilder.)
                     (.name sample)
                     (.alleles [alt-allele])
                     .make)]
              1 [(maybe-fix-vc (first gs) alt-allele)]
              (map :genotype gs)))]
    (if (:prep-sv-genotype config)
      (let [new-gs (ref-vc-genotype (:genotypes orig)
                                    (first (:alt-alleles orig)))]
        (-> orig
            (assoc :vc
              (-> (VariantContextBuilder. (:vc orig))
                  (.genotypes new-gs)
                  .make))
            (assoc :genotypes (map from-genotype new-gs))))
      orig)))

(defn- ordered-vc-iter
  "Provide VariantContexts ordered by chromosome and normalized."
  [rdr vcf-decoder sample config]
  (->> rdr
       line-seq
       (#(if (:prep-sort-pos config) (sort-by-position %) %))
       (remove nochange-alt?)
       (map vcf-decoder)
       (remove #(no-call-genotype? % config))
       (map (partial normalize-sv-genotype config sample))
       (map (partial fix-vc sample config))
       (map :vc)))

(defn- fix-vcf-line
  "Provide fixes to VCF input lines that do not require VariantContext parsing.
  Fixes:
    - INFO lines with empty attributes (starting with ';'), found in
      Complete Genomics VCF files
    - Chromosome renaming."
  [line chr-map config]
  (letfn [(empty-attribute-info [info]
            (if (.startsWith info ";")
              (subs info 1)
              info))
          (fix-info [xs]
            (assoc xs 7 (empty-attribute-info (nth xs 7))))
          (fix-chrom [new xs]
            (assoc xs 0 new))]
    (let [parts (string/split line #"\t")
          cur-chrom (or (get chr-map (first parts)) (first parts))]
      {:chrom cur-chrom
       :line (->> parts
                  (fix-chrom cur-chrom)
                  fix-info
                  (string/join "\t"))})))

(defn- vcf-by-chrom
  "Split input VCF into separate files by chromosome, returning a map of file names."
  [vcf-file ref-file orig-ref-file tmp-dir config]
  (letfn [(ref-chr-files [ref-file]
            (into (ordered-map)
                  (map (fn [x] [x (str (fs/file tmp-dir (str "prep" x ".vcf")))])
                       (chrs-from-fasta-file ref-file))))
          (write-by-chrom [ref-wrtrs chr-map line]
            (let [line-info (fix-vcf-line line chr-map config)]
              (if-let [wtr (get ref-wrtrs (:chrom line-info))]
                (.write wtr (str (:line line-info) "\n"))
                (throw (Exception. (format "Could not find remapping of chromosome %s in reference: %s"
                                           (:chrom line-info) (keys ref-wrtrs)))))))]
    (let [ref-chrs (ref-chr-files ref-file)
          ref-wrtrs (zipmap (keys ref-chrs) (map writer (vals ref-chrs)))
          chr-map (chr-name-remap (:prep-org config) ref-file orig-ref-file)]
      (with-open [rdr (reader vcf-file)]
        (->> rdr
             line-seq
             (drop-while #(.startsWith % "#"))
             (map (partial write-by-chrom ref-wrtrs chr-map))
             doall)
        (doseq [x (vals ref-wrtrs)]
          (.close x)))
      ref-chrs)))

;; ## Top level functionality to manage inputs and writing.

(defn- update-header
  "Update header information, removing contig and adding sample names."
  [sample config]
  (letfn [(clean-metadata [header]
            (apply ordered-set (remove #(= "contig" (.getKey %)) (.getMetaDataInInputOrder header))))]
    (fn [_ header]
      (let [cur-samples (.getGenotypeSamples header)
            new-samples (if (and (:fix-sample-header config)
                                 (not (nil? sample))
                                 (< (count cur-samples) 2))
                          (ordered-set sample)
                          cur-samples)]
        (VCFHeader. (clean-metadata header) new-samples)))))

(defn fix-vcf-sample
  "Update a VCF file with one item to have the given sample name."
  [in-file sample ref]
  (let [out-file (itx/add-file-part in-file "samplefix")]
    (when (itx/needs-run? out-file)
      (with-open [vcf-iter (get-vcf-iterator in-file ref)]
        (write-vcf-w-template in-file {:out out-file}
                              (map #(:vc (fix-vc sample {} %)) (parse-vcf vcf-iter))
                              ref :header-update-fn (update-header sample {}))))
    out-file))

(defn- write-prepped-vcf
  "Write VCF file with correctly ordered and cleaned variants."
  [vcf-file out-info ref-file orig-ref-file sample config]
  (itx/with-temp-dir [tmp-dir (fs/parent (:out out-info))]
    (let [reader-by-chr (into (ordered-map) (map (fn [[k v]] [k (reader v)])
                                                 (vcf-by-chrom vcf-file ref-file orig-ref-file
                                                               tmp-dir config)))]
      (with-open [vcf-reader (AsciiLineReader. (input-stream vcf-file))]
        (let [vcf-decoder (get-vcf-line-parser vcf-reader)]
          (write-vcf-w-template vcf-file out-info
                                (flatten
                                 (for [rdr (vals reader-by-chr)]
                                   (ordered-vc-iter rdr vcf-decoder sample config)))
                                ref-file
                                :header-update-fn (update-header sample config))))
      (doseq [x (vals reader-by-chr)]
        (.close x)))))

(defn prep-vcf
  "Prepare VCF for comparison by normalizing high level attributes
  Assumes by position sorting of variants in the input VCF. Chromosomes do
  not require a specific order, but positions internal to a chromosome do.
  Currently configured for human preparation."
  [in-vcf-file ref-file sample & {:keys [out-dir out-fname config orig-ref-file]
                                  :or {config {}}}]
  (let [config (merge-with #(if (nil? %1) %2 %1) config
                           {:prep-org :GRCh37 :prep-allele-count 2
                            :prep-sort-pos false :prep-sv-genotype false
                            :fix-sample-header false
                            :remove-refcalls true})
        base-name (if (nil? out-fname) (itx/remove-zip-ext in-vcf-file) out-fname)
        out-file (itx/add-file-part base-name "prep" out-dir)]
    (when (itx/needs-run? out-file)
      (write-prepped-vcf in-vcf-file {:out out-file}
                         ref-file orig-ref-file
                         sample config))
    out-file))

(defn pick-best-ref
  "Choose a reference genome for a variant file from set of choices."
  [vcf refs]
  (letfn [(get-vcf-contig [fname]
            (with-open [rdr (reader fname)]
              (->> (line-seq rdr)
                   (drop-while #(.startsWith % "#"))
                   first
                   (#(string/split % #"\t"))
                   first)))
          (has-contig? [contig ref-file]
            (contains?
             (set (keys (get-seq-name-map ref-file)))
             contig))]
    (let [test-contig (get-vcf-contig vcf)]
      (first (filter (partial has-contig? test-contig) refs)))))

;; ## Remove problem characters
;; Handle cleanup for VCF files before feeding to any verifying
;; parser.

(defn- ref-base-getter
  "Given a VCF line, retrieve the reference base. Second optional
   function allows you adjust which base is retrieved to get previous
   bases for indels missing reference padding."
  [ref-file]
  (let [chr-map (prep-rename-map :GRCh37 ref-file)
        get-ref-chrom (fn [chrom]
                        (get chr-map chrom chrom))]
    (fn [xs adjust-fn extra-bases]
      (let [i (adjust-fn (Integer/parseInt (second xs)))]
        (string/upper-case
         (str
          (or (extract-sequence ref-file (get-ref-chrom (first xs)) i (+ i extra-bases))
              (extract-sequence ref-file (first xs) i (+ i extra-bases))
              "N")))))))

(defn- maybe-add-indel-pad-base
  "Check reference and alt alleles for lack of a padding base on indels.
  The VCF spec requires this and GATK will parse incorrectly when a variant
  lacks a shared padding base for indels."
  [ref-file ref-base-get xs]
  (letfn [(get-ref-alts [xs]
            {:ref (nth xs 3) :alts (string/split (nth xs 4) #",")})
          (is-alt-sv? [alt]
            (or (.startsWith alt "<")
                (.contains alt "[")
                (.contains alt "]")))
          (indel? [a]
            (some #(and (not (is-alt-sv? %))
                        (not= (count (:ref a)) (count %)))
                  (:alts a)))
          (is-5pad-n? [a]
            (every? #(and (.startsWith (:ref a) "N") (.startsWith % "N")) (:alts a)))
          (fix-5pad-n [xs a]
            (let [base (ref-base-get xs identity 0)]
              (-> xs
                  (assoc 3 (str base (subs (:ref a) 1)))
                  (assoc 4 (string/join ","
                                        (map #(str base (subs % 1)) (:alts a)))))))
          (no-pad? [a]
            (some #(and (not (is-alt-sv? %))
                        (not= (first (:ref a)) (first %)))
                  (:alts a)))
          (fix-nopad [xs a]
            (let [base (ref-base-get xs dec 0)]
              (-> xs
                  (assoc 1 (dec (Integer/parseInt (second xs))))
                  (assoc 3 (str base (:ref a)))
                  (assoc 4 (string/join ","
                                        (map #(str base %) (:alts a)))))))]
    (if (empty? xs) []
        (let [alleles (get-ref-alts xs)]
          (cond
           (and (indel? alleles) (is-5pad-n? alleles)) (fix-5pad-n xs alleles)
           (and (indel? alleles) (no-pad? alleles)) (fix-nopad xs alleles)
           :else xs)))))

(defn- remove-bad-ref
  "Remove calls where the reference base does not match expected reference allele."
  [ref-base-get xs]
  (letfn [(is-bad-ref? [xs]
            (let [check-bases #{\A \C \G \T}
                  vc-ref (nth xs 3)
                  real-ref (ref-base-get xs identity (dec (count vc-ref)))]
              (and (not (nil? real-ref))
                   (every? check-bases (string/upper-case vc-ref))
                   (every? check-bases (string/upper-case real-ref))
                   (not= (string/upper-case vc-ref) (string/upper-case real-ref)))))]
    (cond
     (empty? xs) []
     (is-bad-ref? xs) []
     :else xs)))

(defn- clean-header-line
  "Fixes problematic information in header lines:
  - Handle renaming of sample names to expected.
  - Handles Illumina special case of SNPs with MAXGT and POLY calls.
    Uses the MAXGT calls which make no prior assumptions about polymorphism.
  - Adds sample names for reads without samples.
  - Problem Number specifications with Number=-1
  - Multiple different uses of BKPT flag in Illumina indel and SV files."
  [line sample]
  (letfn [(fix-bad-alt-header [x]
            (str "##ALT=<ID" (string/replace-first x "##ALT=Type" "") ">"))
          (rename-samples [xs want]
            (let [idx (ffirst (filter (fn [[i x]] (.startsWith x want))
                                      (map-indexed vector xs)))]
              (cond
               idx (assoc (vec xs) idx want)
               (contains? #{0 1} (count xs)) [want]
               (.contains (first xs) "_MAXGT") (cons want (rest xs))
               :else xs)))
          (fix-sample-names [x]
            (let [parts (string/split x #"\t")]
              (cond
               (and sample
                    (> (count parts) 8)) (let [[stay-parts samples] (split-at 9 parts)
                                               fix-samples (if (contains? (set samples) sample)
                                                             samples
                                                             (rename-samples samples sample))]
                                           (string/join "\t" (concat stay-parts fix-samples)))
               (= (count parts) 8) (string/join "\t" (concat parts ["FORMAT" sample]))
               :else x)))
          (clean-header [x]
            (cond
             (.startsWith x "##ALT=Type=") (fix-bad-alt-header x)
             (.startsWith x "##FORMAT=<ID=GL,Number=.,Type=String") ""
             (.startsWith x "#CHROM") (fix-sample-names x)
             :else x))
          (fix-bad-info-number [x]
            (string/replace-first x "Number=-1" "Number=."))
          (fix-bad-bkpt [x]
            (string/replace-first x "ID=BKPT,Number=0,Type=Flag" "ID=BKPT,Number=.,Type=String"))]
    (-> line
        fix-bad-info-number
        fix-bad-bkpt
        clean-header)))

(defn- add-missing-genotypes
  "Add genotypes to samples with no-calls. This handles assumed or non-genotyped
   variations like those that come from many structural variant callers."
  [call parts]
  (if (= (count parts) 8)
    (vec (concat parts ["GT" (string/join "/" (repeat (get call :prep-allele-count 2) "1"))]))
    parts))

(defn- remove-problem-alts
  "Remove lines containing alternative alleles that are duplicated, missing or match reference."
  [xs]
  (let [ref (nth xs 3)
        alts (string/split (nth xs 4) #",")]
    (letfn [(has-duplicate-alts? [alts]
              (not= (count alts) (count (set alts))))]
      (cond
       (empty? xs) []
       (some #(= ref %) alts) []
       (some #(= "." %) alts) []
       (has-duplicate-alts? alts) []
       :else xs))))

(defn- remove-non-ascii
  "Remove lines containing non-ascii characters, normally due to corrupted files."
  [xs]
  (letfn [(ascii? [^Character c]
            (CharUtils/isAsciiPrintable c))
          (all-ascii? [s] (every? ascii? s))]
    (if (all-ascii? (string/join " " xs)) xs [])))

(defn- remove-incorrect-qual
  "Remove lines with impossible negative quality values. QUAL is a phred/log scaled score."
  [xs]
  (letfn [(neg-qual? [xs]
            (let [qual (try
                         (Float/parseFloat (nth xs 5))
                         (catch Exception e 1.0))]
              (< qual 0.0)))]
    (cond
     (empty? xs) []
     (neg-qual? xs) []
     :else xs)))

(defn clean-problem-vcf
  "Clean VCF file which GATK parsers cannot handle due to illegal characters.
  Fixes:
    - Gap characters (-) found in REF or ALT indels.
    - Fixes indels without reference padding or N padding.
    - Removes spaces in INFO fields."
  [in-vcf-file ref-file sample call & {:keys [out-dir]}]
  (let [get-ref-base (ref-base-getter ref-file)
        out-file (itx/add-file-part in-vcf-file "preclean" out-dir)]
    (letfn [(remove-gap [n xs]
              (assoc xs n
                     (-> (nth xs n)
                         (string/replace "-" "")
                         (string/replace "." ""))))
            (fix-info-spaces [xs]
              (assoc xs 7
                     (string/replace (nth xs 7) " " "_")))
            (clean-line [line]
              (if (.startsWith line "#")
                (clean-header-line line sample)
                (->> (string/split line #"\t")
                     (remove-gap 3)
                     (remove-gap 4)
                     fix-info-spaces
                     (add-missing-genotypes call)
                     remove-problem-alts
                     remove-non-ascii
                     remove-incorrect-qual
                     (remove-bad-ref get-ref-base)
                     (maybe-add-indel-pad-base ref-file get-ref-base)
                     (string/join "\t"))))]
      (when (itx/needs-run? out-file)
        (itx/with-tx-file [tx-out-file out-file]
          (with-open [rdr (reader in-vcf-file)
                      wtr (writer tx-out-file)]
            (doall
             (map #(.write wtr (str % "\n"))
                  (remove empty? (map clean-line (line-seq rdr))))))))
      out-file)))
