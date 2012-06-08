(ns bcbio.variation.test.utils
  (:use [midje.sweet]
        [bcbio.variation.normalize :only [prep-vcf]]
        [bcbio.variation.utils.cgmetrics]
        [bcbio.variation.utils.summarize :only [vcf-to-table-config]])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]))

(background
 (around :facts
         (let [data-dir (str (fs/file "." "test" "data"))
               config-dir (str (fs/file "." "config"))
               ref (str (fs/file data-dir "GRCh37.fa"))
               vcf (str (fs/file data-dir "gatk-calls.vcf"))
               cg-vcf (str (fs/file data-dir "cg-normalize.vcf"))
               cg-var (str (fs/file data-dir "cg-masterVar-NA19239.tsv"))
               out-cg-var (itx/add-file-part cg-vcf "prep-cgmetrics")
               out-sum-var (str (itx/file-root vcf) "-variantsum.csv")]
           (doseq [x [out-cg-var out-sum-var]]
             (itx/remove-path x))
           ?form)))

(facts "Add Complete Genomics metrics to VCF file."
  (let [ready-vcf (prep-vcf cg-vcf ref "NA12939" :config {:prep-sort-pos true})]
    (add-cgmetrics ready-vcf cg-var ref) => out-cg-var))

(facts "Summarize variant files, providing tabular output of metrics"
  (let [config-file (str (fs/file config-dir "vcf-summarize.yaml"))]
    (vcf-to-table-config config-file) => [out-sum-var]))
