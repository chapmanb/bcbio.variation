(ns bcbio.variation.test.complex
  "Tests for dealing with more complex variations: structural
  variations and MNPs"
  (:use [midje.sweet]
        [bcbio.variation.complex]
        [bcbio.variation.structural]
        [bcbio.variation.variantcontext])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]))

(let [data-dir (str (fs/file "." "test" "data"))
      ref (str (fs/file data-dir "GRCh37.fa"))
      vcf (str (fs/file data-dir "freebayes-calls-indels.vcf"))
      vc-vcf (str (fs/file data-dir "sv-1000g.vcf"))
      nomnp-out (itx/add-file-part vcf "nomnp")]
  (facts "Deal with multi-nucleotide polymorphisms"
    (normalize-variants vcf ref) => nomnp-out)
  (facts "Parse structural variations"
    (map sv-type (parse-vcf vc-vcf)) => [:DEL :INS :DUP :INS]))
