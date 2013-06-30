(ns bcbio.variation.utils.illumina
  "Automate converting Illumina variant calls into GATK-ready format.
   - Select MAXGT calls from Illumina SNP file (no prior assumption of a variant)
   - Add sample names to SNP and Indel headers.
   - Remove illegal gap characters from indel files.
   - Convert into GRCh37 sorted coordinates.
   - Merge SNP and Indels into single callset."
  (:require [clojure.tools.cli :refer [cli]]
            [me.raynes.fs :as fs]
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
    - orig-ref-file: Original reference file (hg19 for Illumina)
    - out-dir: Output directory to write to.
    - base-tmp-dir: Base temporary directory to work in."
  [base-dir sample-name ref-file orig-ref-file out-dir base-tmp-dir]
  (let [base-dir (fs/expand-home base-dir)
        out-file (str (fs/file out-dir (str sample-name ".vcf")))]
    (when (itx/needs-run? out-file)
      (itx/with-temp-dir [tmp-dir base-tmp-dir]
        (let [call {:name "iprep" :file [(get-illumina-vcf base-dir "SNPs")
                                         (get-illumina-vcf base-dir "Indels")
                                         (get-illumina-vcf base-dir "SVs")]
                    :preclean true :prep true :normalize true
                    :ref orig-ref-file}
              exp {:sample sample-name :ref ref-file}
              out-info (gatk-normalize call exp [] tmp-dir
                                       (fn [_ x] (println x)))]
          (fs/rename (:file out-info) out-file))))
    out-file))

(defn cl-entry [& args]
  (let [[options [base-dir sample-name ref-file orig-ref-file] banner]
        (cli args
             ["-o" "--outdir" "Output directory" :default nil]
             ["-t" "--tmpdir" "Temporary directory (defaults to output director)" :default nil])]
    (when (or (:help options) (some nil? [base-dir sample-name ref-file orig-ref-file]))
      (println "Required arguments:")
      (println "    <base-dir> VCF input file to prepare.")
      (println "    <sample-name> Genome reference file (GRCh37/b37 coordinates)")
      (println "    <ref-file> Genome reference file (GRCh37/b37 coordinates)")
      (println "    <orig-ref-file> Original genome reference file (hg19 coordinates)")
      (println)
      (println banner)
      (System/exit 0))
    (let [out-dir (or (:outdir options) base-dir)
          tmp-dir (or (:tmpdir options) out-dir)
          out-file (prep-illumina-variants base-dir sample-name ref-file orig-ref-file
                                           out-dir tmp-dir)]
      (println out-file)
      (System/exit 0))))
