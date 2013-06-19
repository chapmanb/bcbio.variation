(ns bcbio.variation.test.variantcontext
  (:use [midje.sweet]
        [bcbio.variation.variantcontext])
  (:require [me.raynes.fs :as fs]))

(let [data-dir (str (fs/file "." "test" "data"))
      ref (str (fs/file data-dir "GRCh37.fa"))
      vcf-file (str (fs/file data-dir "gatk-calls.vcf"))]
  (facts "Parsing VCF file to VariantContext"
    (with-open [vcf-iter (get-vcf-iterator vcf-file ref)]
      (let [vc (first (parse-vcf vcf-iter))]
        (:chr vc) => "MT"
        (:start vc) => 73
        (:type vc) => "SNP"
        (:filters vc) => #{}
        (get (:attributes vc) "DP") => "250"
        (-> vc :genotypes count) => 1
        (-> vc :genotypes first :qual) => 99
        (-> vc :genotypes first :type) => "HOM_VAR"
        (-> vc :genotypes first :sample-name) => "Test1"
        (-> vc :genotypes first :attributes) => {"PL" [5820 645 0], "AD" [0 250], "DP" 250
                                                 "GQ" 99 "COV" "0,5"}))))
