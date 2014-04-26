(ns bcbio.variation.utils.comparetwo
  "Compare two variant files within a specified region."
  (:require [clojure.pprint :as pprint]
            [clojure.string :as string]
            [clojure.tools.cli :refer [cli]]
            [clj-yaml.core :as yaml]
            [me.raynes.fs :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.run.fsp :as fsp]
            [bcbio.variation.compare :as compare]
            [bcbio.variation.variantcontext :as gvc]
            [clojure.java.io :as io]))

(defn- prep-compare-config
  "Prepare input configuration for comparing two variant calls in a defined set of regions."
  [eval-vcf refcall-vcf sample ref-file region-file work-dir]
  (let [out-file (str (io/file work-dir "eval-config.yaml"))]
    (->> {:dir {:base work-dir :out "work" :prep "work/prep"}
          :experiments [{:approach "grade" :intervals region-file :ref ref-file :sample sample
                         :calls [{:file refcall-vcf :name "ref" :type "grading-ref" :remove-refcalls true}
                                 {:file eval-vcf :name "eval" :remove-refcalls true}]}]}
         yaml/generate-string
         (spit out-file))
    out-file))

(defn- get-single-sample
  "Retrieve sample from input evaluation file, assuming it is single sample."
  [vcf-file]
  (let [samples (gvc/get-samples vcf-file)]
    (if (= 1 (count samples))
      (first samples)
      (throw (Exception. (str "Expect single sample input VCF. Specify sample with `-s`: " vcf-file))))))

(defn run
  "Run evaluation comparing two inputs"
  [eval-vcf refcall-vcf orig-sample ref-file region-file]
  (let [out-base (fsp/abspath (fsp/file-root (fs/base-name eval-vcf)))
        sample (if (nil? orig-sample) (get-single-sample eval-vcf) orig-sample)
        work-dir (fsp/safe-mkdir (str out-base "-evalwork"))
        config (prep-compare-config (fsp/abspath eval-vcf) (fsp/abspath refcall-vcf)
                                    sample ref-file (fsp/abspath region-file) work-dir)
        conc-out (str out-base "-eval-concordant.vcf")
        disc-out (str out-base "-eval-discordant.vcf")]
    (compare/variant-comparison-from-config config)
    (fs/copy (io/file work-dir "work" (str sample "-eval-ref-concordance.vcf")) conc-out)
    (fs/copy (io/file work-dir "work" (str sample "-eval-ref-discordance.vcf")) disc-out)
    {:concordant conc-out :discordant disc-out}))

(defn cl-entry [& args]
  (let [[options [eval-vcf refcall-vcf ref-file region-file] banner]
        (cli args
             ["-s" "--sample" "Sample to process. If not specified needs single sample input" :default nil])]
    (when (or (:help options) (some nil? [eval-vcf refcall-vcf ref-file region-file]))
      (println "Compare two VCF files in a set of regions defined by a BED file,")
      (println "returning concordants and discordants")
      (println "Required arguments:")
      (println "    <eval-vcf> Evaluation VCF file")
      (println "    <refcall-vcf> Reference call VCF file")
      (println "    <ref-file> Genome reference file (GRCh37/b37 coordinates)")
      (println "    <region-file> BED file of regions to assess in.")
      (println)
      (println banner)
      (System/exit 0))
    (pprint/pprint (run eval-vcf refcall-vcf (:sample options) ref-file region-file))
    (System/exit 0)))
