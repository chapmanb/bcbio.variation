(ns bcbio.variation.test.compare
  (:use [midje.sweet]
        [bcbio.variation.compare])
  (:require [fs.core :as fs]))

(let [data-dir (str (fs/file "." "test" "data"))
      ref (str (fs/file data-dir "hg19.fa"))
      vcf1 (str (fs/file data-dir "gatk-calls.vcf"))
      vcf2 (str (fs/file data-dir "freebayes-calls.vcf"))
      combo-out (add-file-part vcf1 "combine")
      compare-out (str (file-root vcf1) ".eval")]
  (against-background [(before :facts (vec (map #(if (fs/exists? %)
                                                   (fs/delete %))
                                                [combo-out compare-out])))]
    (facts "Variant comparison with GATK"
      (combine-variants vcf1 vcf2 ref) => combo-out
      (variant-comparison vcf1 vcf2 ref) => compare-out)))

(facts "Manipulating file paths"
  (add-file-part "test.txt" "add") => "test-add.txt"
  (add-file-part "/full/test.txt" "new") => "/full/test-new.txt"
  (file-root "/full/test.txt") => "/full/test")
