;; Helper functions to retrieve information from GATK VariantContext
;; objects, which represent variant data stored in VCF files.
(ns bcbio.variation.variantcontext
  (:import [org.broad.tribble.index IndexFactory]
           [org.broad.tribble.source BasicFeatureSource]
           [org.broadinstitute.sting.utils.codecs.vcf VCFCodec])
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
                   (-> vc .getGenotypes .values vec))
   :vc vc})

(defn- vcf-iterator [iter]
  "Lazy iterator over variant contexts in the VCF source iterator."
  (lazy-seq
   (if (.hasNext iter)
     (cons (from-vc (.next iter)) (vcf-iterator iter)))))

(defn parse-vcf [in-file]
  "Lazy iterator of VariantContext information from VCF file."
  (let [idx (IndexFactory/createIndex (file in-file) (VCFCodec.))
        source (BasicFeatureSource. (.getAbsolutePath (file in-file)) idx (VCFCodec.))]
    (vcf-iterator (.iterator source))))
