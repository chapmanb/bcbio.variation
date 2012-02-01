(ns bcbio.variation.test.compare
  (:import [org.broadinstitute.sting.utils.exceptions UserException$BadInput])
  (:use [midje.sweet]
        [bcbio.run.itx]
        [bcbio.variation.annotation]
        [bcbio.variation.callable]
        [bcbio.variation.combine]
        [bcbio.variation.compare]
        [bcbio.variation.filter]
        [bcbio.variation.phasing]
        [bcbio.variation.stats]
        [bcbio.variation.report])
  (:require [fs.core :as fs]))

(let [data-dir (str (fs/file "." "test" "data"))
      ref (str (fs/file data-dir "hg19.fa"))
      intervals (str (fs/file data-dir "target-regions.bed"))
      vcf1 (str (fs/file data-dir "gatk-calls.vcf"))
      vcf2 (str (fs/file data-dir "freebayes-calls.vcf"))
      align-bam (str (fs/file data-dir "aligned-reads.bam"))
      sample "Test1"
      callable-out (format "%s-callable.bed" (file-root align-bam))
      annotated-out (add-file-part vcf2 "annotated")
      combo-out (add-file-part vcf1 "combine")
      compare-out (str (file-root vcf1) ".eval")
      filter-out (add-file-part vcf1 "filter")
      nofilter-out (add-file-part filter-out "nofilter")
      combine-out [(add-file-part vcf1 "fullcombine-wrefs")
                   (add-file-part vcf2 "fullcombine-wrefs")]
      match-out {:concordant (add-file-part combo-out "concordant")
                 :discordant (add-file-part combo-out "discordant")}
      select-out (doall (map #(str (fs/file data-dir (format "%s-%s.vcf" sample %)))
                             ["gatk-freebayes-concordance"
                              "gatk-freebayes-discordance"
                              "freebayes-gatk-discordance"]))]
  (against-background [(before :facts (vec (map #(if (fs/exists? %)
                                                   (fs/delete %))
                                                (concat
                                                 [combo-out compare-out callable-out
                                                  annotated-out filter-out nofilter-out]
                                                 combine-out
                                                 (vals match-out)
                                                 select-out))))]
    (facts "Variant comparison and assessment with GATK"
      (select-by-concordance sample {:name "gatk" :file vcf1}
                             {:name "freebayes" :file vcf2} ref
                             :interval-file intervals) => select-out
      (combine-variants [vcf1 vcf2] ref) => combo-out
      (variant-comparison sample vcf1 vcf2 ref
                          :interval-file intervals) => compare-out
      (-> (concordance-report-metrics sample compare-out)
          first :percent_non_reference_sensitivity) => "88.89"
      (identify-callable align-bam ref) => callable-out
      (let [is-callable? (callable-checker align-bam ref)]
        (is-callable? "chrM" 16 17) => true
        (is-callable? "chrM" 252 252) => false
        (is-callable? "chrM" 5100 5200) => false
        (is-callable? "chrM" 16 15) => false)
      (add-variant-annotations vcf2 align-bam ref) => annotated-out)
    (facts "Create merged VCF files for comparison"
      (create-merged [vcf1 vcf2] [align-bam align-bam] [true true] ref) => combine-out)
    (facts "Filter variant calls avoiding false positives."
      (variant-filter vcf1 ["QD < 2.0" "MQ < 40.0"] ref) => filter-out
      (remove-cur-filters filter-out ref) => nofilter-out
      (split-variants-by-match vcf1 vcf2 ref) => match-out
      (variant-recalibration-filter vcf1 [{:file (:concordant match-out)
                                           :name "concordant"
                                           :prior 10.0}]
                                    ref) => (throws UserException$BadInput
                                                    (contains "annotations with zero variance")))))

(let [data-dir (str (fs/file "." "test" "data"))
      vcf1 (str (fs/file data-dir "gatk-calls.vcf"))]
  (facts "Accumulate statistics associated with variations."
    (map :metric (vcf-stats vcf1)) => ["AC" "AF" "AN" "BaseQRankSum" "DP" "Dels" "FS"
                                       "HRun" "HaplotypeScore" "MQ" "MQ0" "MQRankSum"
                                       "QD" "QUAL" "ReadPosRankSum"]
    (first (vcf-stats vcf1)) => {:max 2.0, :pct75 2.0, :median 2.0, :pct25 2.0, :min 2.0,
                                 :count 10, :metric "AC"}
    (write-summary-table (vcf-stats vcf1)) => nil))

(let [data-dir (str (fs/file "." "test" "data"))
      pvcf (str (fs/file data-dir "phasing-calls.vcf"))
      ref-vcf (str (fs/file data-dir "phasing-reference.vcf"))]
  (facts "Handle haplotype phasing specified in VCF output files."
    (let [haps (parse-phased-haplotypes pvcf)]
      (count haps) => 4
      (count (first haps)) => 5
      (-> haps first first :start) => 10
      (count (second haps)) => 1
      (-> haps (nth 2) first :start) => 16))
  (facts "Compare phased calls to haploid reference genotypes."
    (let [cmps (score-phased-calls pvcf ref-vcf)]
      (map :variant-type (first cmps)) => [:snp :snp :indel :snp :snp]
      (map :comparison (last cmps)) => [:concordant :phasing-error :concordant :discordant]
      (map :nomatch-het-alt (first cmps)) => [true false true false true])))

(facts "Determine the highest count of items in a list"
  (highest-count []) => nil
  (highest-count ["a" "a"]) => "a"
  (highest-count ["a" "b" "b"]) => "b"
  (highest-count ["a" "a" "b" "b"]) => "a"
  (highest-count ["b" "b" "a" "a"]) => "a")
