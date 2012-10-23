(ns bcbio.variation.multisample
  "Compare multiple sample input files, allowing flexible configuration
  of concordance/discordance logic for comparison two sets of calls."
  (:use [clojure.java.io])
  (:require [bcbio.run.itx :as itx]
            [bcbio.variation.variantcontext :as gvc]))

(defn multiple-samples?
  "Check if the input VCF file has multiple genotyped samples."
  [in-file & {:keys [sample]}]
  (let [samples (-> in-file gvc/get-vcf-header .getGenotypeSamples)]
    (or (> (count samples) 1)
        (and (not (nil? sample))
             (not (contains? (set samples) sample))))))

(defn get-out-basename
  "Retrieve basename for output display, handling multiple sample inputs."
  [exp call in-files]
  (let [sample-name (or (:sample exp)
                        (-> in-files first gvc/get-vcf-header .getGenotypeSamples first
                            (str "multi")))]
    (format "%s-%s" sample-name (:name call))))

(defn- write-concordance-output
  [vc-info c1 c2 exp config]
  (let [out-dir (get-in config [:dir :out])
        base-out (str (file out-dir (format "%s-%s.vcf"
                                            (get-out-basename exp c1)
                                            (:name c2))))
        out-files (into {:concordant (itx/add-file-part base-out "concordant")}
                        (map (fn [c]
                               [(keyword (str (:name c) "-discordant"))
                                (itx/add-file-part base-out (str (:name c) "-discordant"))])
                             [c1 c2]))]))

(defn compare-two-vcf-multisample
  "Compare two multisample variant input files.
   TODO: intervals"
  [c1 c2 exp config]
  (let [c2-seen (atom #{})]
    (with-open [c2-get (gvc/get-vcf-retriever (:ref exp) (:file c2))
                c1-iter (gvc/get-vcf-iterator (:file c1) (:ref exp))]
      ))
  (throw (Exception. "Not implemented")))