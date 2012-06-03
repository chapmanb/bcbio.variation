(ns bcbio.variation.test.utils
  (:use [midje.sweet]
        [bcbio.variation.normalize :only [prep-vcf]]
        [bcbio.variation.utils.cgmetrics])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]))

(background
 (around :facts
         (let [data-dir (str (fs/file "." "test" "data"))
               ref (str (fs/file data-dir "GRCh37.fa"))
               cg-vcf (str (fs/file data-dir "cg-normalize.vcf"))
               cg-var (str (fs/file data-dir "cg-masterVar-NA19239.tsv"))
               out-cg-var (itx/add-file-part cg-vcf "prep-cgmetrics")]
           (doseq [x [out-cg-var]]
             (itx/remove-path x))
           ?form)))

(facts "Add Complete Genomics metrics to VCF file."
  (let [ready-vcf (prep-vcf cg-vcf ref "NA12939" :config {:prep-sort-pos true})]
    (add-cgmetrics ready-vcf cg-var ref) => out-cg-var))
