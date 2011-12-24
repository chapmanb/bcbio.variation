(ns bcbio.variation.test.compare
  (:use [midje.sweet]
        [bcbio.variation.compare]))

(facts "Variant comparison with GATK"
  (combine-variants "test1.vcf" "test2.vcf" "") => "test1-combine.vcf")

(facts "Manipulating file paths"
  (add-file-part "test.txt" "add") => "test-add.txt"
  (add-file-part "/full/test.txt" "new") => "/full/test-new.txt"
  (file-root "/full/test.txt") => "/full/test")
