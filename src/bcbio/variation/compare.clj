;; Generate comparisons between two sets of variant calls
;; Utilizes GATK walkers to generate detailed and summary statistics
;; about two sets of calls
;; - Identify non-callable regions with CallableLociWalker
;; - Combine variants from two samples
;; - Use VariantEval to calculate overall concordance statistics
;; - Provide output for concordant and discordant regions for
;;   detailed investigation

(ns bcbio.variation.compare
  (:import [org.broadinstitute.sting.gatk CommandLineGATK])
  (:use [bcbio.variation.variantcontext :only [parse-vcf write-vcf-w-template]]
        [bcbio.variation.stats :only [vcf-stats print-summary-table]]
        [clojure.math.combinatorics :only [combinations]])
  (:require [fs.core :as fs]
            [clj-yaml.core :as yaml]))

;; Utility functions for processing file names

(defn file-root [fname]
  "Retrieve file name without extension: /path/to/fname.txt -> /path/to/fname"
  (let [i (.lastIndexOf fname ".")]
    (if (pos? i)
      (subs fname 0 i)
      fname)))

(defn add-file-part [fname part]
  "Add file extender: base.txt -> base-part.txt"
  (format "%s-%s%s" (file-root fname) part (fs/extension fname)))

(defn- run-gatk [program args]
  (let [std-args ["-T" program "--phone_home" "NO_ET"]]
    (CommandLineGATK/start (CommandLineGATK.)
                           (into-array (concat std-args args)))))

;; GATK walker based variance assessment

(defn combine-variants [vcf1 vcf2 ref]
  "Combine two variant files with GATK CombineVariants."
  (letfn [(unique-name [f]
            (-> f fs/base-name file-root))]
    (let [out-file (add-file-part vcf1 "combine")
          args ["-R" ref
                (str "--variant:" (unique-name vcf1)) vcf1
                (str "--variant:" (unique-name vcf2)) vcf2
                "-o" out-file
                "--genotypemergeoption" "UNIQUIFY"]]
      (if-not (fs/exists? out-file)
        (run-gatk "CombineVariants" args))
      out-file)))

(defn variant-comparison [vcf1 vcf2 ref & {:keys [out-base]
                                           :or {out-base nil}}]
  "Compare two variant files with GenotypeConcordance in VariantEval"
  (let [out-file (str (file-root (if (nil? out-base) vcf1 out-base)) ".eval")
        args ["-R" ref
              "--out" out-file
              "--eval" vcf1
              "--comp" vcf2
              "--evalModule" "GenotypeConcordance"
              "--stratificationModule" "Sample"]]
    (if-not (fs/exists? out-file)
      (run-gatk "VariantEval" args))
    out-file))

(defn select-by-concordance [sample call1 call2 ref]
  "Variant comparison producing 3 files: concordant and both directions discordant"
  (let [base-dir (fs/parent (:file call1))]
    (doall
     (for [[c1 c2 cmp-type] [[call1 call2 "concordance"]
                             [call1 call2 "discordance"]
                             [call2 call1 "discordance"]]]
       (let [out-file (str (fs/file base-dir (format "%s-%s-%s-%s.vcf"
                                                     sample (:name c1) (:name c2) cmp-type)))
             args ["-R" ref
                   "--variant" (:file c1)
                   (str "--" cmp-type) (:file c2)
                   "--out" out-file]]
         (if-not (fs/exists? out-file)
           (run-gatk "SelectVariants" args))
         out-file)))))

;; Custom parsing and combinations using GATK VariantContexts

(defn- vc-by-match-category [in-file]
  "Lazy stream of VariantContexts categorized by concordant/discordant matching."
  (letfn [(genotype-alleles [g]
            (vec (map #(.toString %) (:alleles g))))
          (is-concordant? [vc]
            (= (-> (map genotype-alleles (:genotypes vc))
                   set
                   count)
               1))]
    (for [vc (parse-vcf in-file)]
      [(if (is-concordant? vc) :concordant :discordant)
       vc])))

(defn split-variants-by-match [vcf1 vcf2 ref]
  "Provide concordant and discordant variants for two variant files."
  (let [combo-file (combine-variants vcf1 vcf2 ref)
        out-map {:concordant (add-file-part combo-file "concordant")
                 :discordant (add-file-part combo-file "discordant")}]
    (if-not (fs/exists? (:concordant out-map))
      (write-vcf-w-template combo-file out-map (vc-by-match-category combo-file)
                            ref))
    out-map))

(defn -main [config-file]
  (let [config (-> config-file slurp yaml/parse-string)]
    (doseq [exp (:experiments config)]
      (doseq [[c1 c2] (combinations (:calls exp) 2)]
        (let [c-files (select-by-concordance (:sample exp) c1 c2 (:ref exp))]
          (variant-comparison (:file c1) (:file c2) (:ref exp)
                              :out-base (first c-files))
          (doseq [f c-files]
            (println (fs/base-name f))
            (print-summary-table (vcf-stats f))))))))
