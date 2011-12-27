(ns bcbio.variation.test.compare
  (:use [midje.sweet]
        [bcbio.variation.compare]
        [bcbio.variation.stats])
  (:require [fs.core :as fs]))

(let [data-dir (str (fs/file "." "test" "data"))
      ref (str (fs/file data-dir "hg19.fa"))
      vcf1 (str (fs/file data-dir "gatk-calls.vcf"))
      vcf2 (str (fs/file data-dir "freebayes-calls.vcf"))
      combo-out (add-file-part vcf1 "combine")
      compare-out (str (file-root vcf1) ".eval")
      match-out {:concordant (add-file-part combo-out "concordant")
                 :discordant (add-file-part combo-out "discordant")}
      select-out (doall (map #(str (fs/file data-dir %)) ["gatk-freebayes-concordance.vcf"
                                                          "gatk-freebayes-discordance.vcf"
                                                          "freebayes-gatk-discordance.vcf"]))]
  (against-background [(before :facts (vec (map #(if (fs/exists? %)
                                                   (fs/delete %))
                                                (concat
                                                 [combo-out compare-out]
                                                 (vals match-out)
                                                 select-out))))]
    (facts "Variant comparison with GATK"
      (select-by-concordance "gatk" vcf1 "freebayes" vcf2 ref) => select-out
      (combine-variants vcf1 vcf2 ref) => combo-out
      (variant-comparison vcf1 vcf2 ref) => compare-out
      (split-variants-by-match vcf1 vcf2 ref) => match-out)))

(let [data-dir (str (fs/file "." "test" "data"))
      vcf1 (str (fs/file data-dir "gatk-calls.vcf"))]
  (facts "Accumulate statistics associated with variations."
    (first (vcf-stats vcf1)) => {:max 2.0, :pct75 2.0, :median 2.0, :pct25 2.0, :min 2.0,
                                 :count 10, :metric "AC"}
    (print-summary-table (vcf-stats vcf1)) => nil))

(facts "Manipulating file paths"
  (add-file-part "test.txt" "add") => "test-add.txt"
  (add-file-part "/full/test.txt" "new") => "/full/test-new.txt"
  (file-root "/full/test.txt") => "/full/test")
