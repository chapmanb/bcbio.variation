(ns bcbio.align.interval
  "Convert BED interval contig names between compatible assemblies
  Handles Human hg19 to GRCh37 naming conversions."
  (:use [clojure.java.io]
        [bcbio.align.ref :only [get-seq-name-map]]
        [bcbio.variation.normalize :only [prep-rename-map]])
  (:require [clojure.string :as string]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]))

(defn- update-contig-name
  "Remap the input contig into GRCh37 contig name."
  [name-map line]
  (let [parts (string/split line #"\t")
        remap-contig (get name-map (first parts))]
    (when remap-contig
      (cons remap-contig (rest parts)))))

(defn rename-bed
  "Rename BED coordinates to match supplied reference file"
  [bed-file ref-file & {:keys [out-dir]}]
  (let [out-file (fsp/add-file-part bed-file "remap" out-dir)]
    (when (itx/needs-run? out-file)
      (let [name-map (prep-rename-map :GRCh37 ref-file)]
        (with-open [rdr (reader bed-file)
                    wtr (writer out-file)]
          (doall
           (map #(.write wtr (str (string/join "\t" %) "\n"))
                (->> (line-seq rdr)
                     (map (partial update-contig-name name-map))
                     (remove nil?)))))))
    out-file))
