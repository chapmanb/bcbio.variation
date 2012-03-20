(ns bcbio.align.ref
  "Deal with reference sequences for alignment and variant calling."
  (:import [org.broadinstitute.sting.gatk.datasources.reference ReferenceDataSource]
           [net.sf.picard.reference ReferenceSequenceFileFactory])
  (:use [clojure.java.io]))

(defn get-seq-dict
  "Retrieve Picard sequence dictionary from FASTA reference file."
  [ref-file]
  (ReferenceDataSource. (file ref-file))
  (-> ref-file
      file
      ReferenceSequenceFileFactory/getReferenceSequenceFile
      .getSequenceDictionary))
