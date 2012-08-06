(ns bcbio.variation.multiple
  "Handle useful comparisons from multiple variation calling approaches.
  High level API to consolidate pairwise variant comparisons."
  (:use [clojure.set :only [union]]
        [ordered.map :only [ordered-map]]
        [bcbio.variation.annotation :only [add-variant-annotations]]
        [bcbio.variation.callable :only [get-callable-checker is-callable? has-callers?]]
        [bcbio.variation.combine :only [combine-variants]]
        [bcbio.variation.metrics :only [nonref-passes-filter?]]
        [bcbio.variation.variantcontext :only [parse-vcf get-vcf-retriever
                                               variants-in-region
                                               get-vcf-iterator write-vcf-w-template]])
  (:require [clojure.string :as string]
            [fs.core :as fs]
            [bcbio.run.itx :as itx]
            [bcbio.run.broad :as broad]))

;; ## Utility functions

(defn remove-mod-name [x & {:keys [mods] :or {mods ["recal"]}}]
  "Removes modification names from an approach name."
  (reduce (fn [final mod]
            (string/replace final (str "-" mod) ""))
          x mods))

(defn prep-cmp-name-lookup
  "Lookup map of comparisons by method names.
   - ignore: a list of method names to ignore when creating the lookup map.
   - remove-mods?: Flag to remove naming modifications. This
                   will replace original comparisons with recalibrated."
  [cmps & {:keys [ignore remove-mods?] :or {ignore #{}}}]
  (reduce (fn [m x]
            (let [names (map #(let [n (get-in x [% :name])]
                                (if-not remove-mods? n
                                        (remove-mod-name n :mods [(get-in x [% :mod])])))
                             [:c1 :c2])]
              (if (some #(contains? ignore %) names) m
                  (assoc m names x))))
          (ordered-map)
          cmps))

(defn- not-target? [target-name xs]
  (not (contains? (set (map remove-mod-name xs)) target-name)))

;; ## Prepare multi-overlap sets

(defn get-vc-set-calls
  "Retrieve all called items from a variant context 'set' attribute."
  [vc calls]
  (when-let [set-val (get-in vc [:attributes "set"])]
    (if (= set-val "Intersection")
      (set (map :name calls))
      (->> (string/split set-val #"-")
           (remove #(.startsWith % "filter"))
           (map #(string/split % #"AND"))
           (apply concat)
           set))))

(defmulti select-variant-by-set
  "Select samples based on name of a 'set' from CombineVariants."
  (fn [_ _ set-name] (keyword set-name)))

(defmethod select-variant-by-set :Intersection
  [vcf-in ref set-name]
  (let [file-info {:out-vcf (itx/add-file-part vcf-in set-name nil)}
        args ["-R" ref
              "-o" :out-vcf
              "--variant" vcf-in
              "-select" (format "set == '%s'" set-name)]]
    (broad/run-gatk "SelectVariants" args file-info {:out [:out-vcf]})
    (:out-vcf file-info)))

(defmethod select-variant-by-set :default
  ^{:doc "Select non-intersection names, handling GATK special cases like intersection
          and filtered."}
  [vcf-in ref set-name]
  (letfn [(in-set? [vc]
            (contains? (get-vc-set-calls vc [{:name set-name}])
                       set-name))]
    (let [out-file (itx/add-file-part vcf-in set-name nil)]
      (when (itx/needs-run? out-file)
        (with-open [in-iter (get-vcf-iterator vcf-in ref)]
          (write-vcf-w-template vcf-in {:out out-file}
                                (map :vc (filter in-set? (parse-vcf in-iter)))
                                ref)))
      out-file)))

(defn- gen-all-concordant
  "Create VCF of the intersection of all concordant calls."
  [cmps-by-name ref out-dir config & {:keys [do-include? base-ext]
                                  :or {base-ext "multiall"}}]
  (let [concordant-map (reduce (fn [m [k v]]
                                 (if (or (nil? do-include?) (do-include? k))
                                   (assoc m (get-in v [:c-files :concordant]) (string/join "AND" k))
                                   m))
                               (ordered-map) cmps-by-name)
        union-vcf (combine-variants (keys concordant-map) ref :merge-type :full :out-dir out-dir
                                    :name-map concordant-map :base-ext base-ext)]
    {:union union-vcf
     :intersection (select-variant-by-set union-vcf ref "Intersection")}))

(defmulti gen-target-fps
  "Generate false positives, dispatching differently when the target is from recalling."
  (fn [_ _ _ _ is-recalled? _ _]
    (if is-recalled? :recall :default)))

(defmethod gen-target-fps :recall
  ^{:doc "False positive generation for combine call sets resulting from recalling.
          Report calls found in only a single original set"}
  [target-cmps target-name _ target-overlap-vcf is-recalled? ref out-dir]
  (let [calls (vec (set (mapcat (juxt :c1 :c2) (vals target-cmps))))
        out-file (itx/add-file-part target-overlap-vcf "singles" out-dir)]
    (letfn [(is-single-fp? [vc]
              (= 1 (count (disj (get-vc-set-calls vc calls) target-name))))]
      (when (itx/needs-run? out-file)
        (with-open [in-iter (get-vcf-iterator target-overlap-vcf ref)]
          (write-vcf-w-template target-overlap-vcf {:out out-file}
                                (map :vc (filter is-single-fp? (parse-vcf in-iter)))
                                ref))))
    out-file))

(defmethod gen-target-fps :default
  ^{:doc "False positive generation for single call sets: report discordant variants
          callable in other samples."}
  [target-cmps target-name other-conc-vcf _ _ ref out-dir]
  (letfn [(check-shared [fetch any-callable]
            (fn [x]
              (and (nonref-passes-filter? x)
                   (if (has-callers? any-callable)
                     (is-callable? any-callable (:chr x) (:start x) (:end x))
                     (not (empty? (variants-in-region fetch (:chr x)
                                                      (:start x) (:end x))))))))
          (get-shared-discordant [xs fetch any-callable]
            (let [pass-and-shared? (check-shared fetch any-callable)]
              (map :vc (filter pass-and-shared? xs))))]
    (let [disc-vcfs (remove nil? (map (fn [v]
                                        (get-in v [:c-files
                                                   (keyword (format "%s-discordant" target-name))]))
                                      (vals target-cmps)))
          disc-vcf (-> (combine-variants disc-vcfs ref :merge-type :full :out-dir out-dir
                                         :base-ext (format "dis%s" target-name))
                       (select-variant-by-set ref "Intersection"))
          out-file (itx/add-file-part disc-vcf "shared")
          align-bams (->> (vals target-cmps)
                          (map (juxt :c1 :c2))
                          flatten
                          (map :align)
                          (remove nil?))]
      (with-open [disc-iter (get-vcf-iterator disc-vcf ref)
                  other-retriever (get-vcf-retriever ref other-conc-vcf)
                  call-source (get-callable-checker align-bams ref
                                                    :out-dir (str (fs/parent out-dir)))]
        (when (itx/needs-run? out-file)
          (write-vcf-w-template disc-vcf {:out out-file}
                                (get-shared-discordant (parse-vcf disc-iter)
                                                       other-retriever call-source)
                                ref)))
      out-file)))

(defn- gen-target-problems
  "Create files of false negatives and positives from target-name."
  [target-name target-call cmps-by-name true-p-vcf target-overlap-vcf ref out-dir config]
  (let [notarget-concordant (gen-all-concordant cmps-by-name ref out-dir config
                                                :do-include? (partial not-target? target-name)
                                                :base-ext (format "multino%s" target-name))]
    {:false-negatives
     (-> (combine-variants [true-p-vcf (:intersection notarget-concordant)]
                           ref :merge-type :full :out-dir out-dir
                           :name-map {true-p-vcf "truep"
                                      (:intersection notarget-concordant) target-name}
                           :base-ext (format "multiall-no%s" target-name))
         (select-variant-by-set ref target-name)
         (add-variant-annotations (:align target-call) ref target-call :out-dir out-dir))
     :false-positives (gen-target-fps (remove #(not-target? target-name (first %))
                                              cmps-by-name)
                                      target-name (:union notarget-concordant)
                                      target-overlap-vcf (:recall target-call)
                                      ref out-dir)}))

(defn multiple-overlap-analysis
  "Provide high level concordance overlap comparisons for multiple call approaches.
  Organizes relative to the given target name generating:
   - VCF of calls concordant in all methods: intersection of all concordant calls.
     These are true positives.
   - VCF of calls discordant in the target method, but concordant in the remainder:
     the intersection of all concordant pairs not including target-name minus the
     overall intersection of concordants. These are false negatives.
   - VCF of non-ref calls discordant in the target method and called in any of the other
     methods. We restrict to shared calls to avoid penalizing unique calls.
     These are false positives."
  [cmps config target-name & {:keys [dirname ignore] :or {dirname "multiple"
                                                          ignore #{}}}]
  (let [cmps-by-name (prep-cmp-name-lookup (if (map? cmps) (vals cmps) cmps)
                                           :ignore (union ignore #{"all" "validate"}))
        out-dir (str (fs/file (get-in config [:dir :prep] (get-in config [:dir :out]))
                              dirname))
        ref (-> cmps-by-name vals first :exp :ref)
        target-call (->> cmps-by-name
                         (remove #(not-target? target-name (first %)))
                         first
                         second
                         ((juxt :c1 :c2))
                         (filter #(= (remove-mod-name (:name %)) target-name))
                         first)]
    (when-not (fs/exists? out-dir)
      (fs/mkdirs out-dir))
    (let [all-overlap (gen-all-concordant cmps-by-name ref out-dir config)
          true-p-vcf (add-variant-annotations (:intersection all-overlap) (:align target-call)
                                              ref target-call :out-dir out-dir)
          target-overlaps (-> all-overlap
                              :union
                              (select-variant-by-set ref target-name)
                              (add-variant-annotations (:align target-call) ref
                                                       target-call :out-dir out-dir))
          target-problems (gen-target-problems target-name target-call cmps-by-name
                                               true-p-vcf target-overlaps ref out-dir config)]
      (ordered-map :true-positives true-p-vcf
                   :false-negatives (:false-negatives target-problems)
                   :false-positives (:false-positives target-problems)
                   :target-overlaps target-overlaps))))

(defn pipeline-compare-multiple
  "Perform high level pipeline comparison of a target with multiple experiments."
  [cmps finalizer exp config]
  (let [analysis (multiple-overlap-analysis cmps config (:target finalizer)
                                            :ignore (set (get finalizer :ignore #{})))]
    {:c-files analysis
     :c1 {:name (:target finalizer)}
     :c2 {:name "all"}
     :exp exp :dir (config :dir)}))
