(ns bcbio.variation.utils.sanger
  "Parse and organize results from Sanger validation into VCF"
  (:import [org.broadinstitute.sting.utils.variantcontext Allele
            VariantContextBuilder GenotypeBuilder GenotypesContext]
           [org.broadinstitute.sting.utils.codecs.vcf
            VCFHeader VCFInfoHeaderLine VCFHeaderLineCount VCFHeaderLineType
            VCFFormatHeaderLine]
           [org.broadinstitute.sting.utils.variantcontext.writer
            VariantContextWriterFactory])
  (:require [clojure.string :as string]
            [clojure.java.io :as io]
            [incanter.excel :as excel]
            [incanter.core :as icore]
            [fs.core :as fs]
            [bcbio.align.ref :refer [get-seq-dict]]
            [bcbio.run.itx :as itx]))

(defn- row->varinfo
  "Map a row to variant information"
  [key ref alt group for rev]
  (let [val (cond
             (and (.startsWith for "Pass")
                  (.startsWith rev "Pass")
                  (= for rev))
             (-> for
                 (string/replace "Pass_" "")
                 (string/replace "," ""))
             (and (= for "Fail") (= rev "Fail")) ref
             :else nil)]
    (when val
      (let [[chrom start] (string/split key #"_")
            alts (string/split alt #",")]
        {:chrom chrom
         :start (Integer/parseInt start)
         :ref ref
         :alts (if (= ref val) alts (set (cons val alts)))
         :group group
         :val val}))))

(defn- read-sanger-xls
  "Return confirmed calls from Sanger calls in input XLS file."
  [in-file]
  (->> (excel/read-xls (str in-file))
       (icore/$map row->varinfo ["key" "vcf ref" "vcf allele" "set"
                                 "forward validated" "reverse validated"])
       (remove nil?)))

(defn- get-header [sample-name]
  (VCFHeader. #{(VCFInfoHeaderLine. "valgroup" 1 VCFHeaderLineType/String
                                    "Validation group")
                (VCFFormatHeaderLine. "GT" 1 VCFHeaderLineType/String
                                      "Genotype")}
              #{sample-name}))

(defn- chrom-coord
  "Retrieve chromosome and real position"
  [ref-dict]
  (let [name->index (into {} (map-indexed (fn [i x] [(.getSequenceName x) i])
                                          (.getSequences ref-dict)))]
    (fn [x]
      [(name->index (:chrom x)) (:start x)])))

(defn- val->vc
  "Convert a dictionary of validation information into a VariantContext."
  [sample-name x]
  (println x)
  (-> (VariantContextBuilder. (:chrom x) (:chrom x)
                              (inc (:start x)) (+ (count (:ref x)) (:start x))
                              (cons (Allele/create (:ref x) true)
                                    (map #(Allele/create % false) (:alts x))))
      (.attributes {"valgroup" (:group x)})
      (.genotypes (GenotypesContext/create
                   (java.util.ArrayList.
                    [(GenotypeBuilder/create sample-name
                                             [(Allele/create (:val x) (= (:ref x) (:val x)))])])))
      .make))

(defn sanger->vcf
  "Read directory of Sanger Illumina validation results in Excel format.
   Convert into VCF with validated and reference variants."
  [sanger-dir sample-name ref-file]
  (let [seq-dict (get-seq-dict ref-file)
        out-file (str (fs/file sanger-dir (str sample-name "-sanger-validate.vcf")))]
    (when true ;(itx/needs-run? out-file)
      (with-open [writer (VariantContextWriterFactory/create (io/file out-file)
                                                             seq-dict)]
        (.writeHeader writer (get-header sample-name))
        (doseq [vc (->> (mapcat read-sanger-xls (fs/glob (fs/file sanger-dir "*.xls*")))
                        (sort-by (chrom-coord seq-dict))
                        (map (partial val->vc sample-name)))]
          (.add writer vc))))))