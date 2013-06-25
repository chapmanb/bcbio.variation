(ns bcbio.variation.test.utils
  (:use [midje.sweet]
        [bcbio.variation.annotate.entropy]
        [bcbio.variation.normalize :only [prep-vcf]]
        [bcbio.variation.utils.cgmetrics]
        [bcbio.variation.utils.summarize :only [vcf-to-table-config]])
  (:require [me.raynes.fs :as fs]
            [bcbio.run.itx :as itx]))

(background
 (around :facts
         (let [data-dir (str (fs/file "." "test" "data"))
               config-dir (str (fs/file "." "config"))
               ref (str (fs/file data-dir "GRCh37.fa"))
               alt-ref (str (fs/file data-dir "hg19.fa"))
               vcf (str (fs/file data-dir "gatk-calls.vcf"))
               cg-vcf (str (fs/file data-dir "cg-normalize.vcf"))
               cg-vcf-prep (itx/add-file-part cg-vcf "prep")
               cg-var (str (fs/file data-dir "cg-masterVar-NA19239.tsv"))
               out-cg-var (itx/add-file-part cg-vcf "prep-cgmetrics")
               out-sum-var (str (itx/file-root vcf) "-variantsum.csv")]
           (doseq [x [out-cg-var out-sum-var cg-vcf-prep]]
             (itx/remove-path x))
           ?form)))

(facts "Add Complete Genomics metrics to VCF file."
  (let [ready-vcf (prep-vcf cg-vcf ref "NA12939" :config {:prep-sort-pos true}
                            :orig-ref-file nil)]
    (add-cgmetrics ready-vcf cg-var ref) => out-cg-var))

(facts "Summarize variant files, providing tabular output of metrics"
  (let [config-file (str (fs/file config-dir "vcf-summarize.yaml"))]
    (vcf-to-table-config config-file) => [out-sum-var]))

(facts "Calculate sequence entropy, identifying low-complexity regions."
  (seq-entropy "TTTTTTTTTTTTTTTTTTTT") => 0.0
  (seq-entropy "TAGTAGTAGTAGTAGTAGTA") => (roughly 1.5810)
  (seq-entropy "ACACACACATATATATATAT") => (roughly 1.9784)
  (seq-entropy "ATATGTACACACACACACACATATATATATATTTATATAT") => (roughly 2.3988)
  (seq-entropy "TCTTTTTCATCTCTGAAAACAGTGAGAAAATCCTATGCAT") => (roughly 3.5576))
