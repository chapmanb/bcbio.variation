(ns bcbio.variation.custom.nist
  "Explore variant calling from fosmid data against NIST whole genome datasets"
  (:use [bcbio.variation.filter.attr :only [prep-vc-attr-retriever]])
  (:require [clojure.string :as string]
            [criterium.stats :as stats]
            [bcbio.variation.variantcontext :as gvc]))

(defn- get-nist-filter
  "Retrieve filter information associated with NIST variant calls"
  [retriever vc]
  (if-let [nist-vc (first (gvc/variants-in-region retriever vc))]
    (if (empty? (:filters nist-vc))
      :ref-call
      (first (:filters nist-vc)))
    :no-call))

(defn- get-gms-score
  [attr-get vc]
  (let [attr "gms_illumina"]
    (-> (attr-get [attr] vc)
        (get attr))))

(defn- collect-stats-for-discordant
  "Provide high level statistics on discordant calls."
  [fosmid-file nist-file ref-file]
  (let [attr-get (prep-vc-attr-retriever fosmid-file ref-file)]
    (with-open [vrn-iter (gvc/get-vcf-iterator fosmid-file ref-file)
                nist-retriever (gvc/get-vcf-retriever ref-file nist-file)]
      (reduce (fn [coll vc]
                {:filters (let [filt (get-nist-filter nist-retriever vc)]
                            ;; (when (= :ref-call filt)
                            ;;   (println ((juxt :chr :start) vc) (get-gms-score attr-get vc)))
                            (assoc (:filters coll) filt (inc (get-in coll [:filters filt] 0))))
                 :gms (cons (get-gms-score attr-get vc) (:gms coll))})
              {:filters {} :gms []}
              (filter #(= "SNP" (:type %)) (gvc/parse-vcf vrn-iter))))))

(defn- split-filter-name
  "Split NIST filter names into individual components to summarize."
  [x]
  (-> x
      (string/replace "filtered" "")
      (string/split #"Tranche")
      (#(remove empty? %))))

(defn- counts-to-individual-filters
  [xs]
  (letfn [(split-by-filter [[k v]]
            (when-not (keyword? k)
              (partition 2 (interleave (split-filter-name k) (repeat v)))))]
    (reduce (fn [coll [k v]]
              (assoc coll k (+ v (get coll k 0))))
            {}
            (mapcat split-by-filter xs))))

(defn summarize-discordants
  "Summarize discordant calls between fosmid and NIST calls."
  [fosmid-file nist-file ref-file]
  (let [stats (collect-stats-for-discordant fosmid-file nist-file ref-file)
        ready-gms (filter #(< % 100.0) (:gms stats))
        special-filters [:ref-call :no-call]]
    (doseq [k special-filters]
      (println k (get-in stats [:filters k])))
    (doseq [[k v] (counts-to-individual-filters (:filters stats))]
      (println k v))
    (println "GMS" (count (:gms stats)) (count ready-gms)
             (map #(stats/quantile % ready-gms) [0.0 0.25 0.5 0.75 1.0]))))
