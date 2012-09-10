(ns bcbio.align.interval
  "Convert BED interval contig names between compatible assemblies
  Handles Human hg19 to GRCh37 naming conversions."
  (:use [clojure.java.io]
        [bcbio.align.ref :only [get-seq-name-map]]
        [bcbio.variation.normalize :only [hg19-map]])
  (:require [clojure.string :as string]
            [bcbio.run.itx :as itx]))

(defn- fix-non-version-names
  "Convert any non-versioned names into the representative version in ref-dict."
  [base-map ref-dict]
  (letfn [(find-best-match [x check]
            (first (filter #(.startsWith % x) check)))]
    (reduce (fn [coll [k v]]
              (assoc coll k
                     (if (contains? ref-dict v)
                       v
                       (find-best-match v (keys ref-dict)))))
            {} base-map)))

(defn- add-alt-keys
  "Add alternative key variations: currently handles underscore to dash."
  [base-map modtype]
  {:pre [(= modtype :underscore)]}
  (reduce (fn [coll [k v]]
            (-> coll
                (assoc k v)
                (assoc (string/replace k "_" "-") v)))
          {} base-map))

(defn- prep-name-map
  "Fix GRCh37/hg19 name mappings to handle common problem cases."
  [ref-file]
  (-> hg19-map
      (fix-non-version-names (get-seq-name-map ref-file))
      (add-alt-keys :underscore)))

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
  (let [out-file (itx/add-file-part bed-file "remap" out-dir)]
    (when (itx/needs-run? out-file)
      (let [name-map (prep-name-map ref-file)]
        (with-open [rdr (reader bed-file)
                    wtr (writer out-file)]
          (doall
           (map #(.write wtr (str (string/join "\t" %) "\n"))
                (->> (line-seq rdr)
                     (map (partial update-contig-name name-map))
                     (remove nil?)))))))
    out-file))