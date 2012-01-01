(ns bcbio.variation.test.compare
  (:use [midje.sweet]
        [bcbio.variation.compare]
        [bcbio.variation.stats]
        [bcbio.variation.report]
        [bcbio.variation.callable]
        [bcbio.variation.annotation])
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
      match-out {:concordant (add-file-part combo-out "concordant")
                 :discordant (add-file-part combo-out "discordant")}
      select-out (doall (map #(str (fs/file data-dir (format "%s-%s.vcf" sample %)))
                             ["gatk-freebayes-concordance"
                              "gatk-freebayes-discordance"
                              "freebayes-gatk-discordance"]))]
  (against-background [(before :facts (vec (map #(if (fs/exists? %)
                                                   (fs/delete %))
                                                (concat
                                                 [combo-out compare-out callable-out annotated-out]
                                                 (vals match-out)
                                                 select-out))))]
    (facts "Variant comparison and assessment with GATK"
      (select-by-concordance sample {:name "gatk" :file vcf1}
                             {:name "freebayes" :file vcf2} ref
                             :interval-file intervals) => select-out
      (combine-variants vcf1 vcf2 ref) => combo-out
      (variant-comparison sample vcf1 vcf2 ref
                          :interval-file intervals) => compare-out
      (-> (concordance-report-metrics sample compare-out)
          first :percent_non_reference_sensitivity) => "88.89"
      (split-variants-by-match vcf1 vcf2 ref) => match-out
      (identify-callable align-bam ref) => callable-out
      (add-variant-annotations vcf2 align-bam ref) => annotated-out)))

(let [data-dir (str (fs/file "." "test" "data"))
      vcf1 (str (fs/file data-dir "gatk-calls.vcf"))]
  (facts "Accumulate statistics associated with variations."
    (map :metric (vcf-stats vcf1)) => ["AC" "AF" "AN" "BaseQRankSum" "DP" "Dels" "FS"
                                       "HRun" "HaplotypeScore" "MQ" "MQ0" "MQRankSum"
                                       "QD" "QUAL" "ReadPosRankSum"]
    (first (vcf-stats vcf1)) => {:max 2.0, :pct75 2.0, :median 2.0, :pct25 2.0, :min 2.0,
                                 :count 10, :metric "AC"}
    (write-summary-table (vcf-stats vcf1)) => nil))
