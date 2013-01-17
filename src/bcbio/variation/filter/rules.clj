(ns bcbio.variation.filter.rules
  "Define filtration rules used to help identify true/false positives for variant classification.
   Helps organize the logic of selecting variants."
  (:require [bcbio.variation.filter.attr :as attr]
            [bcbio.variation.metrics :as metrics]
            [bcbio.variation.multiple :as multiple]))

(defn vc-zygosity [vc]
  (if (some #(.startsWith (:type %) "HET") (:genotypes vc)) :het :hom))

(defn- below-support-thresh?
  "Check if a variant context has a low amount of supporting variant calls."
  [vc _ call exp]
  (let [freq (get call :fp-freq 0.25)
        thresh (Math/ceil (* freq (dec (count (:calls exp)))))]
    (-> (multiple/get-vc-set-calls vc (:calls exp))
        (disj (:name call))
        count
        (<= thresh))))

(defn- is-intersection? [vc _ _ _]
  (when-let [set-val (get-in vc [:attributes "set"])]
    (= set-val "Intersection")))

(defn- novel-variant?
  "Is a variant novel, or is it represented in dbSNP?"
  [vc _ _ _]
  (contains? #{nil "."} (:id vc)))

(defn- het-snp? [vc _ _ _]
  (and (= "SNP" (:type vc))
       (= :het (vc-zygosity vc))))

(defn- het-indel? [vc _ _ _]
  (and (not= "SNP" (:type vc))
       (= :het (vc-zygosity vc))))

(defn novel-het-indel? [vc g c e]
  (and (het-indel? vc g c e) (novel-variant? vc g c e)))

(defn- low-call-confidence?
  "Define low confidence calls."
  [vc attr-get _ _]
  (let [attrs (attr-get ["PL" "PLratio"] vc)]
    (when (not (nil? (get attrs "PL")))
      (or (> (get attrs "PL") -7.5)
          (< (get attrs "PLratio" Float/MAX_VALUE) 0.25)))))

(defn- good-pl-support?
  "Identify PL ratios with reasonable support for being a variant."
  [vc attr-get _ _]
  (let [attrs (attr-get ["PLratio"] vc)]
    (when (not-any? nil? (vals attrs))
      (> (get attrs "PLratio") 0.4))))

(defn- low-confidence-novel-het-snp?
  [vc attr-get c e]
  (and (low-call-confidence? vc attr-get c e)
       (novel-variant? vc attr-get c e)
       (het-snp? vc attr-get c e)))

(defn- flex-low-call-confidence?
  "Define calls with a more flexible low confidence call"
  [vc attr-get _ _]
   (let [attrs (attr-get ["PL"] vc)]
    (when (not-any? nil? (vals attrs))
      (> (get attrs "PL") -20.0))))

(defn- low-depth?
  "Calls with low supporting depth"
  [vc attr-get _ _]
  (when (and (het-indel? vc attr-get) (novel-variant? vc))
    (let [attrs (attr-get ["DP"] vc)]
      (when (not-any? nil? (vals attrs))
        (< (get attrs "DP") 25.0)))))

(defn- passes-mapping-quality?
  "Avoid feeding low quality mapping into true/false positives."
  [vc attr-get _ _]
  (let [attrs (attr-get ["MQ"] vc)]
    (when (not-any? nil? (vals attrs))
      (> (get attrs "MQ") 50.0))))

(defn- artifact-allele-balance?
  "Identify skewed allele balances indicative of artifacts.
   This is a signature of problem heterozygote calls from GATK Haplotype caller."
  [vc attr-get _ _]
  (let [attrs (attr-get ["AD"] vc)]
    (when (not-any? nil? (vals attrs))
      (> (get attrs "AD") 0.35))))

(defn- passes-filter?
  [vc _ _ _]
  (metrics/passes-filter? vc))

(def ^{:private true
       :doc "Define keyword mappings to function definitions"}
  rules {:below-call-support below-support-thresh?
         :all-callers is-intersection? 
         :novel novel-variant?
         :het-snp het-snp?
         :het-indel het-indel?
         :novel-het-indel novel-het-indel?
         :low-confidence-novel-het-snp low-confidence-novel-het-snp?
         :low-confidence low-call-confidence?
         :good-pl good-pl-support?
         :flex-low-confidence flex-low-call-confidence?
         :low-depth low-depth?
         :passes-filter passes-filter?
         :high-map-quality passes-mapping-quality?
         :problem-allele-balance artifact-allele-balance?})

(defn vc-checker
  "Identify variants conforming to supplied rules."
  [orig-file call exp]
  (let [attr-get (attr/prep-vc-attr-retriever orig-file (:ref exp))]
    (letfn [(call-rule [vc rulekw]
              ((get rules rulekw) vc attr-get call exp))]
      (fn [vc & {:keys [yes no]}]
        (and (every? (partial call-rule vc) yes)
             (not-any? (partial call-rule vc) no))))))
