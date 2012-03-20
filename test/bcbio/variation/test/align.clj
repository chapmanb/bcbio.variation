(ns bcbio.variation.test.align
  "Test retrieval and analysis of BAM alignment data"
  (:use [midje.sweet]
        [bcbio.align.reorder])
  (:require [fs.core :as fs]))

(let [data-dir (str (fs/file "." "test" "data"))
      ref (str (fs/file data-dir "GRCh37.fa"))
      align-bam (str (fs/file data-dir "aligned-needfix.bam"))]
  (facts "Reorder BAM files to match contigs."
    (chromsome-remap ["22" "MT"] ["22" "MT" "X"]) => []
    (chromsome-remap ["MT" "22"] ["22" "MT" "X"]) =future=> [["22" 9] ["MT" 1]]
    (chromsome-remap ["chrM" "chr22"] ["22" "MT" "X"]) =future=> [["chr22" 9] ["chrM" 1]]
    (reorder-bam align-bam ref {:sample "Test1"} {}) => nil))
