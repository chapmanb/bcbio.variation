(ns bcbio.variation.utils.illumina
  "Automate converting Illumina variant calls into GATK-ready format.
   - Select MAXGT calls from Illumina SNP file (no prior assumption of a variant)
   - Add sample names to SNP and Indel headers.
   - Remove illegal gap characters from indel files.
   - Convert into GRCh37 sorted coordinates.
   - Merge SNP and Indels into single callset."
  (:require [me.raynes.fs :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.variation.combine :refer [gatk-normalize]]))

(defn- get-illumina-vcf
  [base-dir base-name]
  (-> (fs/file base-dir "Variations" (str base-name "*.vcf"))
      str
      fs/glob
      first
      str))

(defn prep-illumina-variants
  "Prepare Illumina variants from a standard directory structure.
    - base-dir: Directory containing Illumina information (will have subdirs like
                Assembly, Consensus and Variations)
    - sample-name: The name to include in updated VCF headers
    - ref-file: Reference file we want to sort to
    - orig-ref-file: Original reference file (hg19 for Illumina)"
  [base-dir sample-name ref-file orig-ref-file]
  (let [base-dir (fs/expand-home base-dir)
        out-file (str (fs/file base-dir "Variations" (str sample-name ".vcf")))]
    (when (itx/needs-run? out-file)
      (itx/with-temp-dir [out-dir base-dir]
        (let [call {:name "iprep" :file [(get-illumina-vcf base-dir "SNPs")
                                         (get-illumina-vcf base-dir "Indels")
                                         (get-illumina-vcf base-dir "SVs")]
                    :preclean true :prep true :normalize true
                    :ref orig-ref-file}
              exp {:sample sample-name :ref ref-file}
              out-info (gatk-normalize call exp [] out-dir
                                       (fn [_ x] (println x)))]
          (fs/rename (:file out-info) out-file))))
    out-file))
