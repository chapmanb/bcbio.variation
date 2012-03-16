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
      nomnp-out (itx/add-file-part vcf "nomnp")]
  (against-background [(before :facts (vec (map itx/remove-path [nomnp-out])))]
    (facts "Deal with multi-nucleotide polymorphisms"
      (normalize-variants vcf ref) => nomnp-out)))

(facts "Parse structural variations"
  (let [data-dir (str (fs/file "." "test" "data"))
        ref (str (fs/file data-dir "GRCh37.fa"))
        vc-vcf (str (fs/file data-dir "sv-1000g.vcf"))]
    (parse-sv-vcf vc-vcf ref) => nil
    (with-open [vcf-source (get-vcf-source vc-vcf ref)]
      (doall (map get-sv-type (parse-vcf vcf-source)))) => (concat (repeat 6 :BND)
                                                                   [nil :DEL :INS :DUP :INS])))
