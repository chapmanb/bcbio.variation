(ns bcbio.variation.test.grade
  "Tests for grading and evaluation of variant calls against a reference set of variations."
  (:use [clojure.java.io]
        [midje.sweet]
        [bcbio.variation.compare :exclude [-main]]
        [bcbio.variation.phasing]
        [bcbio.variation.variantcontext :exclude [-main]])
  (:require [me.raynes.fs :as fs]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.variation.compare :as compare]
            [bcbio.variation.combine :as combine]
            [bcbio.variation.grade :as grade]
            [bcbio.variation.report :as report]
            [bcbio.variation.utils.quickcompare :as qcmp]))

(background
 (around :facts
         (let [data-dir (str (fs/file "." "test" "data"))
               ref-file (str (fs/file data-dir "GRCh37.fa"))]
           (doseq [x []]
             (fsp/remove-path x)
             (when (.endsWith x ".vcf")
               (fsp/remove-path (str x ".idx"))))
           ?form)))

(facts "Compare correct-ploidy sex and mito calls against GIAB reference"
  (let [config-dir (str (fs/file data-dir "giab-comparison"))
        out-dir (str (fs/file config-dir "grading"))
        config-file (str (fs/file config-dir "giab-comparison.yaml"))]
    (fsp/remove-path out-dir)
    (compare/variant-comparison-from-config config-file)))

(facts "Compare diploid phased and haploid callsets."
  (let [base-dir (fs/file data-dir "phased")
        calls {true [{:file (str (fs/file base-dir "NA12878-fosfinal.vcf"))
                      :name "fosfinal"}]
               false [{:file (str (fs/file base-dir "NA12878-illumina.vcf"))
                       :name "illumina"}]}
        exp {:ref ref-file :sample "NA12878"
             :intervals (str (fs/file base-dir "NA12878-cmp-regions.bed"))
             :approach "grade"}
        config {:dir {:out (str (fs/file base-dir "work"))}}]
    (doseq [x (concat [(get-in config [:dir :out])]
                      (fs/glob (str (fs/file base-dir "NA12878-cmp-regions-*")))
                      (fs/glob (str (fs/file base-dir "NA12878-fosfinal-cmp*"))))]
      (fsp/remove-path x))
    (-> (compare-two-vcf-phased calls exp config) :c-files keys) => [:concordant :discordant
                                                                     :discordant-missing :phasing-error]))

(facts "Compare diploid callset against a diploid reference.
        Includes comparison of calls with END tags which can be out of sync following comparison."
  (let [base-dir (fs/file data-dir "digrade")
        c1 {:file (str (fs/file base-dir "NA12878-cmp-r1.vcf"))
            :name "ref" :type "grading-ref"}
        c2 {:file (str (fs/file base-dir "NA12878-cmp-r2.vcf"))
            :name "eval"}
        exp {:ref ref-file :sample "NA12878" :approach "grade"}
        config {:dir {:out (str (fs/file base-dir "work"))}}
        out-file (str (file (get-in config [:dir :out])
                            "NA12878-eval-ref-discordance-annotate.vcf"))]
    (fsp/remove-path (get-in config [:dir :out]))
    (fs/mkdirs (get-in config [:dir :out]))
    (let [cmp (-> [(compare-two-vcf c1 c2 exp config)]
                   (finalize-comparisons exp config)
                   first)]
      (-> cmp :grade-breakdown :discordant :snp :shared :hethom) => 1
      (-> cmp :c-files :eval-discordant) => out-file)))

(facts "Comparisons where inputs have alternative reference sequence: hg19->GRCh37"
  (let [base-dir (fs/file data-dir "digrade")
        c1o {:file (str (fs/file base-dir "NA12878-cmp-r1.vcf"))
             :prep true
             :name "ref" :type "grading-ref"}
        c2o {:file (str (fs/file base-dir "hg19" "NA12878-cmp-r3.vcf"))
             :name "eval"
             :prep true
             :intervals (str (fs/file base-dir "hg19" "NA12878-cmp-r3-intervals.bed"))
             :ref (str (fs/file data-dir "hg19.fa"))}
        exp {:ref ref-file :sample "NA12878" :approach "grade"
             :intervals (str (fs/file base-dir "NA12878-cmp-r1-intervals.bed"))}
        config {:dir {:out (str (fs/file base-dir "work"))}}]
    (fsp/remove-path (get-in config [:dir :out]))
    (fs/mkdirs (get-in config [:dir :out]))
    (let [[c1 c2] (#'compare/prepare-vcf-calls (assoc exp :calls [c1o c2o]) config)]
      (compare-two-vcf c1 c2 exp config))))

(facts "Normalize input VCFs containing END tags"
  (let [base-dir (fs/file data-dir "digrade")
        c1 {:file (str (fs/file base-dir "NA12878-cmp-r1.vcf"))
            :name "ref" :type "grading-ref" :prep true}
        exp {:ref ref-file :sample "NA12878" :approach "grade"}
        out-dir (str (fs/file base-dir "work"))]
    (fsp/remove-path out-dir)
    (fs/mkdirs out-dir)
    (combine/gatk-normalize c1 exp nil out-dir (fn [& xs] (println xs)))))

(facts "Perform quick comparisons between smaller variant files that can fit in memory."
  (let [base-dir (fs/file data-dir "digrade")
        c1 (str (fs/file base-dir "NA12878-cmp-r1.vcf"))
        c2 (str (fs/file base-dir "NA12878-cmp-r2.vcf"))
        dir-out-file (str (fs/file base-dir "NA12878-cmp-r1-cmp.csv"))]
    (fsp/remove-path dir-out-file)
    (qcmp/two-vcfs c1 c2 ref-file) => {:concordant 14 :discordant 3 :sample "NA12878"}
    (qcmp/vcfdir-to-base c1 base-dir ref-file 2) => dir-out-file))
