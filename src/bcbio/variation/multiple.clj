(ns bcbio.variation.multiple
  "Handle useful comparisons from multiple variation calling approaches.
  High level API to consolidate pairwise variant comparisons."
  (:use [ordered.map :only [ordered-map]]
        [bcbio.variation.combine :only [combine-variants]])
  (:require [clojure.string :as string]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]))

(defn prep-cmp-name-lookup
  "Lookup map of comparisons by method names."
  [cmps]
  (reduce (fn [m x] (assoc m [(-> x :c1 :name)
                              (-> x :c2 :name)] x))
          (ordered-map)
          cmps))

(defn- select-variant-by-set
  "Select samples based on name of a 'set' from CombineVariants."
  [vcf-in ref set-name & {:keys [out-dir]}]
  (let [file-info {:out-vcf (itx/add-file-part vcf-in set-name out-dir)}
        args ["-R" ref
              "-o" :out-vcf
              "--variant" vcf-in
              "-select" (format "set == '%s'" set-name)]]
    (broad/run-gatk "SelectVariants" args file-info {:out [:out-vcf]})
    (:out-vcf file-info)))

(defn- gen-all-concordant
  "Create VCF of the intersection of all concordant calls."
  [cmps-by-name out-dir config & {:keys [do-include? base-ext]
                                  :or {base-ext "multiall"}}]
  (let [concordant-map (reduce (fn [m [k v]]
                                 (if (or (nil? do-include?) (do-include? k))
                                   (assoc m (-> v :c-files first) (string/join "-" k))
                                   m))
                               (ordered-map) cmps-by-name)
        ref (-> cmps-by-name vals first :exp :ref)
        out-dir (str (fs/file (get-in config [:dir :prep] (get-in config [:dir :out]))
                              "multiple"))]
    (-> (combine-variants (keys concordant-map) ref :merge-type :full :out-dir out-dir
                          :name-map concordant-map :base-ext base-ext)
        (select-variant-by-set ref "Intersection"))))

(defn- gen-nontarget-concordant
  "Create file of false negatives of concordant samples not called in target-name."
  [target-name cmps-by-name true-p-vcf out-dir config]
  (letfn [(not-target? [x]
            (not (contains? (set x) target-name)))]
    (let [ref (-> cmps-by-name vals first :exp :ref)
          notarget-concordant (gen-all-concordant cmps-by-name out-dir config
                                                  :do-include? not-target?
                                                  :base-ext (format "multino%s" target-name))]
      (-> (combine-variants [true-p-vcf notarget-concordant] ref :merge-type :full :out-dir out-dir
                            :name-map {true-p-vcf "truep" notarget-concordant target-name}
                            :base-ext (format "multiall-no%s" target-name))
          (select-variant-by-set ref target-name)))))

(defn multiple-overlap-analysis
  "Provide high level concordance overlap comparisons for multiple call approaches.
  Organizes relative to the given target name generating:
   - VCF of calls concordant in all methods: intersection of all concordant calls.
     These are true positives.
   - VCF of calls discordant in the target method, but concordant in the remainder:
     the intersection of all concordant pairs not including target-name minus the
     overall intersection of concordants. These are false negatives."
  [cmps config target-name]
  (let [cmps-by-name (prep-cmp-name-lookup cmps)
        out-dir (str (fs/file (get-in config [:dir :prep] (get-in config [:dir :out]))
                              "multiple"))]
    (when-not (fs/exists? out-dir)
      (fs/mkdirs out-dir))
    (let [true-p-vcf (gen-all-concordant cmps-by-name out-dir config)]
      {:true-positives true-p-vcf
       :false-negatives (gen-nontarget-concordant target-name cmps-by-name true-p-vcf
                                                  out-dir config)})))
