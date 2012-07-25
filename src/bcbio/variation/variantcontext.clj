(ns bcbio.variation.variantcontext
  "Helper functions to retrieve information from GATK VariantContext
   objects, which represent variant data stored in VCF files."
  (:import [org.broad.tribble.index IndexFactory]
           [org.broad.tribble AbstractFeatureReader]
           [org.broad.tribble.readers AsciiLineReader]
           [org.broadinstitute.sting.utils.codecs.vcf VCFCodec]
           [org.broadinstitute.sting.utils.variantcontext.writer VariantContextWriterFactory
            Options]
           [org.broadinstitute.sting.gatk.refdata.tracks RMDTrackBuilder]
           [org.broadinstitute.sting.gatk.arguments ValidationExclusion$TYPE]
           [java.util EnumSet])
  (:use [clojure.java.io]
        [lazymap.core :only [lazy-hash-map]]
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
   :attributes (into {} (.getAttributes g))
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

(defprotocol VcfRetrievable
  "Provide a retriever of variants from zero to many inputs."
  (has-variants? [this space start end ref alt])
  (variants-in-region [this space start end]))

(defrecord VcfRetriever [sources]
  VcfRetrievable
  (has-variants? [this space start end ref alt]
    (some #(and (= start (:start %))
                (= end (:end %))
                (= ref (:ref-allele %))
                (= alt (:alt-alleles %)))
          (variants-in-region this space start end)))
  (variants-in-region [_ space start end]
    (mapcat #(map from-vc (iterator-seq (.query % space start end)))
            sources))
  java.io.Closeable
  (close [_]
    (doseq [x sources]
      (.close x))))

(defn get-vcf-retriever
  "Indexed VCF file retrieval for zero to multiple files with clean handle closing."
  [ref & vcf-files]
  (VcfRetriever. (->> vcf-files
                      (remove nil?)
                      (map #(get-vcf-source % ref)))))

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
          (.add (get writer-map fkey) item))
        (doseq [x (vals writer-map)]
          (.close x))))))

(defn write-vcf-from-filter
  "Write VCF file from input using a filter function."
  [vcf ref out-part passes?]
  (with-open [source (get-vcf-source vcf ref)]
    (write-vcf-w-template vcf {:out (itx/add-file-part vcf out-part)}
                          (map :vc (filter passes? (parse-vcf source)))
                          ref)))

(defn -main [vcf ref approach]
  (with-open [vcf-s (get-vcf-source vcf ref)]
    (letfn [(item-iter []
              (case approach
                "line" (map :vc (line-vcf-parser vcf))
                "gatk" (iterator-seq (.iterator vcf-s))
                "orig" (map :vc (parse-vcf vcf-s))))]
      (write-vcf-w-template vcf {:out "vctest.vcf"} (item-iter) ref)
      ;; (doseq [[i x] (map-indexed vector (item-iter))]
      ;;   (when (= 0 (mod i 10000))
      ;;      (println x)))
      )))
