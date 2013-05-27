(ns bcbio.variation.filter.specific
  "Identify technology or caller specific variants from multiple combined callsets."
  (:import [org.broadinstitute.variant.variantcontext VariantContextBuilder]
           [org.broadinstitute.variant.vcf VCFHeader VCFInfoHeaderLine
            VCFHeaderLineCount VCFHeaderLineType])
  (:use [ordered.set :only [ordered-set]]
        [bcbio.variation.filter.attr :only [get-vc-attrs]]
        [bcbio.variation.filter.trusted :only [variant-set-metadata
                                               get-comparison-fullcombine]])
  (:require [bcbio.variation.variantcontext :as gvc])
  (:require [clojure.string :as string]
            [bcbio.run.itx :as itx]))

(defn- get-specific
  [data kw-want kw-cmp]
  (when (and (= 1 (count (kw-want data)))
             (> (count (kw-cmp data)) 1))
    (first (kw-want data))))

(defn get-x-specific-designation
  "Check if a variant is specific to a caller or method."
  [vc calls]
  (let [data (variant-set-metadata vc calls)]
    (reduce (fn [coll [kw-want kw-cmp]]
              (if-let [x (get-specific data kw-want kw-cmp)]
                (assoc coll kw-want x)
                coll))
            {} [[:caller :technology] [:technology :caller]])))

(defn- add-x-specific
  "Add specificity information to a VariantContext if present."
  [vc calls]
  (letfn [(xspec-to-string [[k v]]
            (str (name k) ":" v))]
    (let [xspec (get-x-specific-designation vc calls)]
      (when (seq xspec)
        (-> (VariantContextBuilder. (:vc vc))
            (.attributes (assoc (:attributes vc)
                           "xspecific" (->> xspec
                                            (map xspec-to-string)
                                            (string/join ","))))
            .make)))))

(defn- add-x-specific-header
  [_ header]
  (let [new #{(VCFInfoHeaderLine.
               "xspecific" VCFHeaderLineCount/UNBOUNDED VCFHeaderLineType/String
               "Identify variant call as specific to a technology or calling method.")}]
    (VCFHeader. (apply ordered-set (concat (.getMetaDataInInputOrder header) new))
                (.getGenotypeSamples header))))

(defn poor-call-support?
  "Simple measure to evaluate call support based on depth and allele balance.
   This identifies poorly supported items, which primarily make up false
   positive examples."
  [vc & {:keys [thresh]
         :or {thresh {:dp 100 :ad 0.05}}}]
  (let [attrs (get-vc-attrs vc [[:format "AD"] [:format "DP"]] {})]
    (and (when-let [dp (get attrs [:format "DP"])]
           (< dp (:dp thresh)))
         (when-let [ad (get attrs [:format "AD"])]
           (> ad (:ad thresh))))))

(defn get-x-specific-variants
  "Filter VCF file generating output only variants specific to a technology or caller."
  [cmps support exp config]
  (when-let [base-vcf (get-comparison-fullcombine cmps support config)]
    (let [out-file (itx/add-file-part base-vcf "xspecific")]
      (when (itx/needs-run? out-file)
        (with-open [base-vcf-iter (gvc/get-vcf-iterator base-vcf (:ref exp))]
          (gvc/write-vcf-w-template base-vcf {:out out-file}
                                    (->> (gvc/parse-vcf base-vcf-iter)
                                         (filter poor-call-support?)
                                         (map #(add-x-specific % (:calls exp)))
                                         (remove nil?))
                                    (:ref exp) :header-update-fn add-x-specific-header)))
      out-file)))
