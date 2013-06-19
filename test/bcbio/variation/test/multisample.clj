(ns bcbio.variation.test.multisample
  "Handle variant comparisons of variant files with multiple samples."
  (:use [midje.sweet]
        [clojure.java.io]
        [bcbio.variation.compare :only [variant-comparison-from-config]]
        [bcbio.variation.multisample])
  (:require [me.raynes.fs :as fs]
            [bcbio.run.itx :as itx]))

(background
 (around :facts
         (let [data-dir (str (fs/file "." "test" "data"))
               config-dir (str (fs/file "." "config"))
               config-file-m (str (fs/file config-dir "multiple-comparison.yaml"))
               multi-outdir (str (fs/file data-dir "multi"))
               ref-file (str (fs/file data-dir "GRCh37.fa"))
               vcf-single (str (fs/file data-dir "gatk-calls.vcf"))
               vcf-single2 (str (fs/file data-dir "freebayes-calls.vcf"))
               vcf-m1 (str (fs/file data-dir "1000genome-multi.vcf"))
               vcf-m2 (str (fs/file data-dir "1000genome-multi-cmp.vcf"))
               compare-out {:concordant (str (fs/file multi-outdir "HG00101multi-c1-c2-concordant.vcf"))
                            :c1-discordant (str (fs/file multi-outdir "HG00101multi-c1-c2-c1-discordant.vcf"))
                            :c2-discordant (str (fs/file multi-outdir "HG00101multi-c1-c2-c2-discordant.vcf"))}]
           (doseq [x (concat [] (vals compare-out))]
             (itx/remove-path x)
             (when (.endsWith x ".vcf")
               (itx/remove-path (str x ".idx"))))
           ?form)))

(facts "Check for multiple samples in a VCF file"
  (multiple-samples? vcf-single) => false
  (multiple-samples? vcf-m1) => true)

(facts "Perform comparisons between multiple sample variant files."
  (-> (variant-comparison-from-config config-file-m) first :c-files) => (contains compare-out))

(facts "Perform flexible comparisons between two variant files."
  (let [c1 {:name "gatk" :file vcf-single}
        c2 {:name "freebayes" :file vcf-single2}
        exp {:sample "Test1" :ref ref-file :params {:compare-approach :approximate}}
        config {:dir {:out (str (fs/file data-dir "flex"))}}]
    (itx/remove-path (get-in config [:dir :out]))
    (fs/mkdir (get-in config [:dir :out]))
    (-> (compare-two-vcf-flexible c1 c2 exp config) keys) => (contains :c-files)))