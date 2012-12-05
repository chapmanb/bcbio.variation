(ns bcbio.variation.variantcontext
  "Helper functions to retrieve information from GATK VariantContext
   objects, which represent variant data stored in VCF files."
  (:import [org.broad.tribble.index IndexFactory]
           [org.broad.tribble AbstractFeatureReader]
           [org.broad.tribble.readers AsciiLineReader]
           [org.broadinstitute.sting.utils.codecs.vcf
            VCFCodec VCFUtils VCFHeader VCFFilterHeaderLine]
           [org.broadinstitute.sting.utils.variantcontext VariantContextBuilder]
           [org.broadinstitute.sting.utils.variantcontext.writer VariantContextWriterFactory
            Options]
           [org.broadinstitute.sting.gatk.refdata.tracks RMDTrackBuilder]
           [org.broadinstitute.sting.gatk.arguments ValidationExclusion$TYPE]
           [org.apache.log4j Logger]
           [java.util EnumSet])
  (:use [clojure.java.io]
        [clojure.set :only [intersection union]]
        [lazymap.core :only [lazy-hash-map]]
        [ordered.set :only [ordered-set]]
        [bcbio.align.ref :only [get-seq-dict]])
  (:require [clojure.string :as string]
            [bcbio.run.itx :as itx]))

;; ## Represent VariantContext objects
;;
;; Provide simple map-based access to important attributes of
;; VariantContexts. There are 3 useful levels of abstraction:
;;
;;  - VariantContext: Details about a variation. This captures a
;;    single line in a VCF file
;;  - Genotype: An individual genotype for a sample, at a variant position.
;;  - Allele: The actual alleles at a genotype.

(defn from-genotype
  "Represent a sample genotype including alleles.
   :genotype stores the original java genotype object for direct access."
  [g]
  (lazy-hash-map
   :sample-name (.getSampleName g)
   :qual (.getPhredScaledQual g)
   :type (-> g .getType .name)
   :phased? (.isPhased g)
   :attributes (merge {"DP" (.getDP g) "AD" (vec (.getAD g))
                       "GQ" (.getGQ g) "PL" (vec (.getPL g))}
                      (into {} (.getExtendedAttributes g)))
   :alleles (vec (.getAlleles g))
   :genotype g))

(defn from-vc
  "Provide a top level map of information from a variant context.
   :vc stores the original java VariantContext object for direct access."
  [vc]
  (lazy-hash-map
   :chr (.getChr vc)
   :start (.getStart vc)
   :end (.getEnd vc)
   :id (when (.hasID vc) (.getID vc))
   :ref-allele (.getReference vc)
   :alt-alleles (vec (.getAlternateAlleles vc))
   :type (-> vc .getType .name)
   :filters (set (.getFilters vc))
   :attributes (into {} (.getAttributes vc))
   :qual (.getPhredScaledQual vc)
   :num-samples (.getNSamples vc)
   :genotypes (map from-genotype
                   (-> vc .getGenotypes .toArray vec))
   :vc vc))

;; ## Parsing VCF files

(defn get-vcf-source
  "Create a Tribble FeatureSource for VCF file.
   Handles indexing and parsing of VCF into VariantContexts.
   We treat gzipped files as tabix indexed VCFs."
  [in-file ref-file & {:keys [ensure-safe codec]}]
  (let [cur-codec (if (nil? codec) (VCFCodec.) codec)]
    (if (.endsWith in-file ".gz")
      (AbstractFeatureReader/getFeatureReader in-file cur-codec false)
      (let [validate (when (false? ensure-safe)
                       ValidationExclusion$TYPE/ALLOW_SEQ_DICT_INCOMPATIBILITY)
            idx (.loadIndex (RMDTrackBuilder. (get-seq-dict ref-file) nil validate)
                            (file in-file) cur-codec)]
        (AbstractFeatureReader/getFeatureReader (.getAbsolutePath (file in-file)) cur-codec idx)))))

(defn get-vcf-iterator
  "Create an iterator over VCF VariantContexts."
  [in-file ref-file & {:keys [ensure-safe codec]}]
  (.iterator (get-vcf-source in-file ref-file :ensure-safe ensure-safe
                             :codec codec)))

(defn variants-in-region
  "Retrieve variants located in potentially multiple variant files"
  ([retriever vc]
     (variants-in-region retriever (:chr vc) (:start vc) (:end vc)))
  ([retriever space start end]
     (letfn [(get-vcs-in-source [[source fname]]
               (with-open [vcf-iter (.query source space start end)]
                 (doall (map #(assoc (from-vc %) :fname fname) (iterator-seq vcf-iter)))))]
       (mapcat get-vcs-in-source (map vector (:sources retriever) (:fnames retriever))))))

(defn has-variants?
  "Look for matching variants present in any of the variant files."
  [retriever space start end ref alt]
  (some #(and (= start (:start %))
              (= end (:end %))
              (= ref (:ref-allele %))
              (seq (intersection (set (:alt-alleles %)) (set alt))))
        (variants-in-region retriever space start end)))

(defrecord VariantRetriever [sources fnames]
  java.io.Closeable
  (close [_]
    (doseq [x sources]
      (.close x))))

(defn get-vcf-retriever
  "Indexed variant file retrieval for zero to multiple files with clean handle closing."
  [ref & vcf-files]
  (let [fnames (remove nil? vcf-files)]
    (VariantRetriever. (map #(get-vcf-source % ref) fnames)
                       fnames)))

(defn parse-vcf
  "Lazy iterator of VariantContext information from VCF file."
  [vcf-source]
  (map from-vc (iterator-seq (.iterator vcf-source))))

(defn get-vcf-line-parser
  "Retrieve parser to do line-by-line parsing of VCF files."
  [vcf-reader]
  (let [codec (VCFCodec.)]
    (.readHeader codec vcf-reader)
    (fn [line]
      (from-vc (.decode codec line)))))

(defn- line-vcf-parser
  [vcf]
  (let [parser (with-open [rdr (AsciiLineReader. (input-stream vcf))]
                 (get-vcf-line-parser rdr))]
    (map parser (drop-while #(.startsWith % "#") (line-seq (reader vcf))))))

(defn get-vcf-header
  "Retrieve header from input VCF file."
  [vcf-file]
  (with-open [vcf-reader (AsciiLineReader. (input-stream vcf-file))]
    (.readHeader (VCFCodec.) vcf-reader)))

;; ## Writing VCF files

(defn merge-headers
  [& merge-files]
  (fn [_ header]
    (VCFHeader. (VCFUtils/smartMergeHeaders (cons header (map get-vcf-header merge-files))
                                            (Logger/getLogger ""))
                (.getGenotypeSamples header))))

(defn write-vcf-w-template
  "Write VCF output files starting with an original input template VCF.
   Handles writing to multiple VCF files simultaneously with the different
   file handles represented as keywords. This allows lazy splitting of VCF files:
   `vc-iter` is a lazy sequence of `(writer-keyword variant-context)`.
   `out-file-map` is a map of writer-keywords to output filenames."
  [tmpl-file out-file-map vc-iter ref & {:keys [header-update-fn]}]
  (letfn [(make-vcf-writer [f ref]
            (VariantContextWriterFactory/create (file f) (get-seq-dict ref)
                                                (EnumSet/of Options/INDEX_ON_THE_FLY
                                                            Options/ALLOW_MISSING_FIELDS_IN_HEADER)))
          (convert-to-output [info]
            [(if (and (coll? info) (= 2 (count info))) (first info) :out)
             (if (coll? info) (last info) info)])]
    (itx/with-tx-files [tx-out-files out-file-map (keys out-file-map) [".idx"]]
      (let [tmpl-header (get-vcf-header tmpl-file)
            writer-map (zipmap (keys tx-out-files)
                               (map #(make-vcf-writer % ref) (vals tx-out-files)))]
        (doseq [[key out-vcf] writer-map]
          (.writeHeader out-vcf (if-not (nil? header-update-fn)
                                  (header-update-fn key tmpl-header)
                                  tmpl-header)))
        (doseq [[fkey item] (map convert-to-output vc-iter)]
          (let [ready-vc (if (contains? item :vc) (:vc item) item)]
            (when-not (nil? ready-vc)
              (.add (get writer-map fkey) ready-vc))))
        (doseq [x (vals writer-map)]
          (.close x))))))

(defn- add-filter-header
  [fname fdesc]
  (fn [_ header]
    (let [new #{(VCFFilterHeaderLine. fname fdesc)}]
      (VCFHeader. (apply ordered-set (concat (.getMetaDataInInputOrder header) new))
                  (.getGenotypeSamples header)))))

(defn- maybe-add-filter
  [fname passes? vc]
  (if (passes? vc)
    (:vc vc)
    (-> (VariantContextBuilder. (:vc vc))
        (.filters (union #{fname} (:filters vc)))
        .make)))

(defn write-vcf-from-filter
  "Write VCF file from input using a filter function."
  [vcf ref out-part fname fdesc passes?]
  (let [out-file (itx/add-file-part vcf out-part)]
    (when (itx/needs-run? out-file)
      (with-open [vcf-iter (get-vcf-iterator vcf ref)]
        (write-vcf-w-template vcf {:out out-file}
                              (map (partial maybe-add-filter fname passes?) (parse-vcf vcf-iter))
                              ref
                              :header-update-fn (add-filter-header fname fdesc))))
    out-file))

(defn select-variants
  "Select variants from an input file with supplied filter."
  [in-file passes? file-out-part ref-file]
  (let [out-file (itx/add-file-part in-file file-out-part)]
    (when (itx/needs-run? out-file)
      (with-open [in-iter (get-vcf-iterator in-file ref-file)]
        (write-vcf-w-template in-file {:out out-file}
                              (map :vc (filter passes? (parse-vcf in-iter)))
                              ref-file)))
    out-file))

(defn -main [vcf ref approach]
  (with-open [vcf-iter (get-vcf-iterator vcf ref)]
    (letfn [(item-iter []
              (case approach
                "line" (map :vc (line-vcf-parser vcf))
                "gatk" (iterator-seq (.iterator vcf-iter))
                "orig" (map :vc (parse-vcf vcf-iter))))]
      (write-vcf-w-template vcf {:out "vctest.vcf"} (item-iter) ref)
      ;; (doseq [[i x] (map-indexed vector (item-iter))]
      ;;   (when (= 0 (mod i 10000))
      ;;      (println x)))
      )))
