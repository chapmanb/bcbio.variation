;; Helper functions to retrieve information from GATK VariantContext
;; objects, which represent variant data stored in VCF files.
(ns bcbio.variation.variantcontext
  (:import [org.broad.tribble.index IndexFactory]
           [org.broad.tribble.source BasicFeatureSource]
           [org.broadinstitute.sting.utils.codecs.vcf VCFCodec StandardVCFWriter]
           [net.sf.picard.reference ReferenceSequenceFileFactory])
  (:use [clojure.java.io]))

(defn- from-genotype [g]
  "Top level map of useful details from a genotype."
  {:sample-name (.getSampleName g)
   :qual (.getPhredScaledQual g)
   :type (-> g .getType .name)
   :attributes (into {} (.getAttributes g))
   :alleles (vec (.getAlleles g))
   :genotype g})

(defn from-vc [vc]
  "Provide a top level map of information from a variant context."
  {:chr (.getChr vc)
   :start (.getStart vc)
   :end (.getEnd vc)
   :ref-allele (.getReference vc)
   :type (-> vc .getType .name)
   :filters (set (.getFilters vc))
   :attributes (into {} (.getAttributes vc))
   :genotypes (map from-genotype
                   (-> vc .getGenotypes .toArray vec))
   :vc vc})

(defn- vcf-iterator [iter]
  "Lazy iterator over variant contexts in the VCF source iterator."
  (lazy-seq
   (if (.hasNext iter)
     (cons (from-vc (.next iter)) (vcf-iterator iter)))))

(defn- vcf-source [in-file]
  (let [idx (IndexFactory/createIndex (file in-file) (VCFCodec.))]
    (BasicFeatureSource. (.getAbsolutePath (file in-file)) idx (VCFCodec.))))

(defn parse-vcf [in-file]
  "Lazy iterator of VariantContext information from VCF file."
  (vcf-iterator (.iterator (vcf-source in-file))))

(defmacro with-open-map [binding-map & body]
  "Emulate with-open with bindings supplied as a map."
  `(try
     ~@body
     (finally
      (vec (map #(.close %) (vals ~binding-map))))))

(defn write-vcf-w-template [tmpl-file out-file-map vc-iter ref]
  "Write VCF output files starting with an original input template VCF.
  - vc-iter is a lazy sequence of (writer-keyword variant-context)
  - out-file-map is a map of writer-keywords to output filenames."
  (letfn [(make-vcf-writer [f ref]
            (StandardVCFWriter. (file f)
                                (.getSequenceDictionary
                                 (ReferenceSequenceFileFactory/getReferenceSequenceFile (file ref)))))]
    (let [tmpl-header (.getHeader (vcf-source tmpl-file))
          writer-map (zipmap (keys out-file-map)
                             (map #(make-vcf-writer % ref) (vals out-file-map)))]
      (with-open-map writer-map
        (doseq [out-vcf (vals writer-map)]
          (.writeHeader out-vcf tmpl-header))
        (doseq [[category vc] vc-iter]
          (.add (get writer-map category) (:vc vc)))))))
