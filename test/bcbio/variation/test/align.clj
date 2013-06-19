(ns bcbio.variation.test.align
  "Test retrieval and analysis of BAM alignment data"
  (:use [midje.sweet]
        [bcbio.align.reorder])
  (:require [me.raynes.fs :as fs]
            [bcbio.run.itx :as itx]))

(let [data-dir (str (fs/file "." "test" "data"))
      ref (str (fs/file data-dir "GRCh37.fa"))
      align-bam (str (fs/file data-dir "aligned-needfix.bam"))
      out-bam (itx/add-file-part align-bam "reorder")]
  (against-background [(before :facts (doall (map itx/remove-path
                                                  [out-bam (str align-bam ".bai")])))]
    (facts "Reorder BAM files to match contigs."
      (get-new-chr-order ["22" "MT"] ["22" "MT" "X"] ref) => nil
      (get-new-chr-order ["MT" "22"] ["22" "MT" "X"] ref) => {:names ["22" "MT"]
                                                              :indexes {1 0 0 1}}
      (get-new-chr-order ["chrM" "chr22"] ["22" "MT" "X"] ref ) => {:names ["chr22" "chrM"]
                                                                    :indexes {1 0 0 1}}
      (reorder-bam align-bam ref {} {:sample "Test1"}) => out-bam)))
