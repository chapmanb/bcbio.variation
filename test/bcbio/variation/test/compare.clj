(ns bcbio.variation.test.compare
  (:import [org.broadinstitute.sting.utils.exceptions UserException$BadInput])
  (:use [midje.sweet]
        [bcbio.run.itx]
        [bcbio.variation.annotation]
        [bcbio.variation.callable]
        [bcbio.variation.combine]
        [bcbio.variation.compare]
        [bcbio.variation.filter]
        [bcbio.variation.normalize]
        [bcbio.variation.phasing]
        [bcbio.variation.metrics]
        [bcbio.variation.multiple]
        [bcbio.variation.report]
        [bcbio.variation.variantcontext])
  (:require [fs.core :as fs]))

(let [data-dir (str (fs/file "." "test" "data"))
      ref (str (fs/file data-dir "GRCh37.fa"))
      intervals (str (fs/file data-dir "target-regions.bed"))
      vcf1 (str (fs/file data-dir "gatk-calls.vcf"))
      vcf2 (str (fs/file data-dir "freebayes-calls.vcf"))
      align-bam (str (fs/file data-dir "aligned-reads.bam"))
      sample "Test1"
      annotated-out (add-file-part vcf2 "annotated")
      combo-out (add-file-part vcf1 "combine")
      compare-out (str (file-root vcf1) ".eval")
      filter-out (add-file-part vcf1 "filter")
      nofilter-out (add-file-part filter-out "nofilter")
      combine-out [(add-file-part vcf1 "fullcombine-wrefs")
                   (add-file-part vcf2 "fullcombine-wrefs")]
      combine-out-xtra [(add-file-part vcf1 "mincombine")
                        (add-file-part vcf1 "fullcombine")
                        (add-file-part vcf2 "fullcombine")]
      match-out {:concordant (add-file-part combo-out "concordant")
                 :discordant (add-file-part combo-out "discordant")}
      select-out (doall (map #(str (fs/file data-dir (format "%s-%s.vcf" sample %)))
                             ["gatk-freebayes-concordance"
                              "gatk-freebayes-discordance"
                              "freebayes-gatk-discordance"]))]
  (against-background [(before :facts (vec (map #(if (fs/exists? %)
                                                   (fs/delete %))
                                                (concat
                                                 [combo-out compare-out 
                                                  annotated-out filter-out nofilter-out]
                                                 combine-out combine-out-xtra
                                                 (vals match-out)
                                                 select-out))))]
    (facts "Variant comparison and assessment with GATK"
      (select-by-concordance sample {:name "gatk" :file vcf1}
                             {:name "freebayes" :file vcf2} ref
                             :interval-file intervals) => select-out
      (combine-variants [vcf1 vcf2] ref) => combo-out
      (calc-variant-eval-metrics sample vcf1 vcf2 ref
                                 :intervals intervals) => compare-out
      (-> (concordance-report-metrics sample compare-out)
          first :percent_non_reference_sensitivity) => "88.89"
      (add-gatk-annotations vcf2 align-bam ref) => annotated-out)
    (facts "Create merged VCF files for comparison"
      (create-merged [vcf1 vcf2] [align-bam align-bam] [true true] ref) => combine-out)
    (facts "Filter variant calls avoiding false positives."
      (variant-filter vcf1 ["QD < 2.0" "MQ < 40.0"] ref) => filter-out
      (remove-cur-filters filter-out ref) => nofilter-out
      (split-variants-by-match vcf1 vcf2 ref) => match-out
      (variant-recal-filter vcf1 [{:file (:concordant match-out)
                                   :name "concordant"
                                   :truth "true"
                                   :prior 10.0}]
                            ["QD" "DP"]
                            ref) => (throws UserException$BadInput
                                            (contains "Error during negative model training")))))

(let [data-dir (str (fs/file "." "test" "data"))
      ref (str (fs/file data-dir "GRCh37.fa"))
      align-bam (str (fs/file data-dir "aligned-reads.bam"))
      out-callable (format "%s-callable.bed" (file-root align-bam))
      out-intervals (add-file-part out-callable "intervals")]
  (against-background [(before :facts (doall (map remove-path
                                                  [out-callable out-intervals])))]
    (facts "Check for callability based on sequencing reads."
      (identify-callable align-bam ref) => out-callable
      (let [[is-callable? call-source] (callable-checker align-bam ref)]
        (with-open [_ call-source]
          (is-callable? "MT" 16 17) => true
          (is-callable? "MT" 252 252) => false
          (is-callable? "MT" 5100 5200) => false
          (is-callable? "MT" 16 15) => false))
      (get-callable-bed align-bam ref) => out-intervals)))

(facts "Accumulate statistics associated with variations."
  (let [data-dir (str (fs/file "." "test" "data"))
        ref (str (fs/file data-dir "GRCh37.fa"))
        vcf1 (str (fs/file data-dir "gatk-calls.vcf"))
        vcf2 (str (fs/file data-dir "freebayes-calls.vcf"))]
    (map :metric (vcf-stats vcf1 ref)) => ["AC" "AF" "AN" "BaseQRankSum" "DP" "Dels" "FS"
                                           "HRun" "HaplotypeScore" "MQ" "MQ0" "MQRankSum"
                                           "QD" "QUAL" "ReadPosRankSum"]
    (first (vcf-stats vcf1 ref)) => {:max 2.0, :pct75 2.0, :median 2.0, :pct25 2.0, :min 2.0,
                                     :count 10, :metric "AC"}
    (write-summary-table (vcf-stats vcf1 ref)) => nil
    (let [metrics (get-vcf-classifier-metrics ref [vcf1 vcf2])
          wnil-metrics (get-vcf-classifier-metrics ref [vcf1 vcf2] :remove-nil-cols false)]
      (count metrics) => 2
      (-> metrics first :cols) => ["AC" "AF" "AN" "DP" "QUAL"]
      (-> metrics first :rows count) => 10
      (-> wnil-metrics first :cols) => ["AC" "AF" "AN" "DP" "QUAL" "ReadPosRankSum"]
      (-> wnil-metrics first :rows count) => 2
      (-> metrics second :rows first) => [2.0 1.0 2.0 938.0 99.0]
      (classify-decision-tree metrics) => ["DP"]
      (classify-decision-tree wnil-metrics) => []
      (merge-classified-metrics [["A" "B" "C"] ["C" "D"]]) => {:top-metrics ["A" "C" "B" "D"]})))

(let [data-dir (str (fs/file "." "test" "data"))
      ref (str (fs/file data-dir "GRCh37.fa"))
      pvcf (str (fs/file data-dir "phasing-contestant.vcf"))
      ref-vcf (str (fs/file data-dir "phasing-reference.vcf"))
      ref2-vcf (str (fs/file data-dir "phasing-reference2.vcf"))]
  (facts "Handle haplotype phasing specified in VCF output files."
    (with-open [pvcf-source (get-vcf-source pvcf ref)]
      (let [haps (parse-phased-haplotypes pvcf-source)]
        (count haps) => 4
        (count (first haps)) => 6
        (-> haps first first :start) => 9
        (count (second haps)) => 1
        (-> haps (nth 2) first :start) => 16)))
  (facts "Compare phased calls to haploid reference genotypes."
    (with-open [ref-vcf-s (get-vcf-source ref-vcf ref)
                pvcf-s (get-vcf-source pvcf ref)]
      (let [cmps (score-phased-calls pvcf-s ref-vcf-s)]
        (map :variant-type (first cmps)) => [:snp :snp :indel :snp :snp]
        (:comparison (ffirst cmps)) => :discordant
        (map :comparison (last cmps)) => [:ref-concordant :phasing-error
                                          :ref-concordant :discordant]
        (map :nomatch-het-alt (first cmps)) => [false true false false true])))
  (facts "Compare two sets of haploid reference calls"
    (with-open [ref-vcf-s (get-vcf-source ref-vcf ref)
                ref2-vcf-s (get-vcf-source ref2-vcf ref)]
      (let [cmps (score-phased-calls ref2-vcf-s ref-vcf-s)]
        (count cmps) => 1
        (count (first cmps)) => 10
        (drop 6 (map :comparison (first cmps))) => [:ref-concordant :concordant
                                                    :ref-concordant :discordant])))
  (facts "Check is a variant file is a haploid reference."
    (is-haploid? pvcf ref) => false
    (is-haploid? ref-vcf ref) => true))

(let [data-dir (str (fs/file "." "test" "data"))
      ref (str (fs/file data-dir "GRCh37.fa"))
      cbed (str (fs/file data-dir "phasing-contestant-regions.bed"))
      rbed (str (fs/file data-dir "phasing-reference-regions.bed"))]
  (facts "Merging and count info for reference and contestant analysis regions."
    (count-comparison-bases rbed cbed ref) => (contains {:compared 12 :total 13})
    (count-comparison-bases rbed nil ref) => (contains {:compared 13 :total 13})))

(facts "Calculate final accuracy score for contestant/reference comparison."
  (calc-accuracy {:total-bases {:compared 10}
                  :discordant {:indel 1 :snp 1}
                  :phasing-error {:indel 1 :snp 1}}) => (roughly 62.50))

(let [data-dir (str (fs/file "." "test" "data"))
      ref (str (fs/file data-dir "GRCh37.fa"))
      vcf (str (fs/file data-dir "cg-normalize.vcf"))
      out-vcf (add-file-part vcf "prep")
      prevcf (str (fs/file data-dir "illumina-needclean.vcf"))
      out-prevcf (add-file-part prevcf "preclean")]
  (against-background [(before :facts (vec (map remove-path [out-vcf out-prevcf
                                                             (str vcf ".idx")])))]
    (facts "Check for multiple samples in a VCF file"
      (multiple-samples? vcf) => false)
    (facts "Normalize variant representation of chromosomes, order, genotypes and samples."
      (prep-vcf vcf ref "Test1" :sort-pos true) => out-vcf)
    (facts "Pre-cleaning of problematic VCF input files"
      (clean-problem-vcf prevcf) => out-prevcf)))

(let [config-file (str (fs/file "." "config" "method-comparison.yaml"))
      config (load-config config-file)
      out-dir (str (fs/file (get-in config [:dir :prep]) "multiple"))
      cmps (variant-comparison-from-config config-file)]
  (letfn [(get-out-files [x ext]
            {:true-positives
             (str (fs/file out-dir (format "Test1-multiall-fullcombine-Intersection%s.vcf" ext)))
             :false-negatives
             (str (fs/file out-dir (format "Test1-multiall-no%s-fullcombine-%s%s.vcf" x x ext)))
             :false-positives
             (str (fs/file out-dir (format "Test1-dis%s-fullcombine-Intersection-shared.vcf" x)))
             :target-overlaps
             (str (fs/file out-dir (format "Test1-multiall-fullcombine-%s%s.vcf" x ext)))})]
    (against-background [(before :facts (vec (map remove-path [out-dir])))]
      (facts "Handle multiple variant approach comparisons."
        (multiple-overlap-analysis cmps config "cg") => (get-out-files "cg" "")
        (multiple-overlap-analysis cmps config "gatk") => (get-out-files "gatk" "-annotated")))))

(facts "Load configuration files, normalizing input."
  (let [config-file (fs/file "." "config" "method-comparison.yaml")
        config (load-config config-file)]
    (get-in config [:dir :out]) => (has-prefix "/")
    (-> config :experiments first :sample) => "Test1"
    (-> config :experiments first :calls first :file) => (has-prefix "/")
    (-> config :experiments first :calls second :filters first) => "HRun > 5.0"))

(facts "Determine the highest count of items in a list"
  (highest-count []) => nil
  (highest-count ["a" "a"]) => "a"
  (highest-count ["a" "b" "b"]) => "b"
  (highest-count ["a" "a" "b" "b"]) => "a"
  (highest-count ["b" "b" "a" "a"]) => "a")
