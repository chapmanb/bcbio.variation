;; Tests for dealing with more complex variations: structural
;; variations and MNPs

(ns bcbio.variation.test.complex
  (:use [midje.sweet]
        [bcbio.variation.complex])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]))

(let [data-dir (str (fs/file "." "test" "data"))
      ref (str (fs/file data-dir "hg19.fa"))
      vcf (str (fs/file data-dir "freebayes-calls-indels.vcf"))
      nomnp-out (itx/add-file-part vcf "nomnp")]
  (facts "Deal with multi-nucleotide polymorphisms"
    (normalize-mnps vcf ref) => nomnp-out))
