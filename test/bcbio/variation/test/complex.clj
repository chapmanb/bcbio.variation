(ns bcbio.variation.test.complex
  "Tests for dealing with more complex variations: structural
  variations and MNPs"
  (:use [midje.sweet]
        [bcbio.variation.complex]
        [bcbio.variation.structural]
        [bcbio.variation.variantcontext])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.itree :as itree]))

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
        vcf (str (fs/file data-dir "sv-1000g.vcf"))
        ill-vcf (str (fs/file data-dir "sv-illumina.vcf"))
        vcf-list (parse-vcf-sv ill-vcf ref)
        vcf-itree (parse-vcf-sv vcf ref :out-format :itree)]
    (-> vcf-list first :start-ci) => 6066065
    (-> vcf-itree (get "22") (itree/iget-range 15883520 15883620) first :end-ci) => 15883626
    (with-open [vcf-source (get-vcf-source vcf ref)
                illvcf-source (get-vcf-source ill-vcf ref)]
      (doall (map get-sv-type (parse-vcf vcf-source))) => (concat (repeat 6 :BND)
                                                                  [nil :DEL :INS :DUP :INS])
      (doall (map get-sv-type (parse-vcf illvcf-source))) => [:DUP :BND :BND :INS :CNV :DEL :INV])))
