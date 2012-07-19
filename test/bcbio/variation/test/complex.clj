(ns bcbio.variation.test.complex
  "Tests for dealing with more complex variations: structural
  variations and MNPs"
  (:use [midje.sweet]
        [bcbio.variation.complex]
        [bcbio.variation.structural]
        [bcbio.variation.variantcontext])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]))

(background
 (around :facts
         (let [data-dir (str (fs/file "." "test" "data"))
               ref (str (fs/file data-dir "GRCh37.fa"))
               target-bed (str (fs/file data-dir "target-regions.bed"))
               sv-vcf1 (str (fs/file data-dir "sv-1000g.vcf"))
               sv-vcf2 (str (fs/file data-dir "sv-illumina.vcf"))
               sv-out {:sv-concordant
                       (str (fs/file data-dir "Test-sv1000g-svIll-svconcordance.vcf"))
                       :sv-sv1000g-discordant
                       (str (fs/file data-dir "Test-sv1000g-svIll-svdiscordance.vcf"))
                       :sv-svIll-discordant
                       (str (fs/file data-dir "Test-svIll-sv1000g-svdiscordance.vcf"))}
               sv-out2
               {:sv-concordant (str (fs/file data-dir "Test-sv1-sv2-svconcordance.vcf"))
                :sv-sv1-discordant (str (fs/file data-dir "Test-sv1-sv2-svdiscordance.vcf"))
                :sv-sv2-discordant (str (fs/file data-dir "Test-sv2-sv1-svdiscordance.vcf"))}
               mnp-vcf (str (fs/file data-dir "freebayes-calls-indels.vcf"))
               cindel-vcf (str (fs/file data-dir "freebayes-calls-complexindels.vcf"))
               cindel-out (itx/add-file-part cindel-vcf "nomnp")
               indel-vcf1 (str (fs/file data-dir "sv-indels-fb.vcf"))
               indel-vcf2 (str (fs/file data-dir "sv-indels-gatk.vcf"))
               indel-out (str (fs/file data-dir "Test-svindfb-svindgatk-svconcordance.vcf"))
               nomnp-out (itx/add-file-part mnp-vcf "nomnp")
               params {:min-indel 100}]
           (doseq [x (concat [nomnp-out indel-out cindel-out] (vals sv-out) (vals sv-out2))]
             (itx/remove-path x))
           ?form)))

(facts "Deal with multi-nucleotide polymorphisms"
  (normalize-variants mnp-vcf ref) => nomnp-out)

(facts "Split complex indels into individual components"
  (normalize-variants cindel-vcf ref) => cindel-out)

(facts "Parse structural variations"
  (let [vcf-list (parse-vcf-sv sv-vcf2 ref)
        vcf-by-region (parse-vcf-sv sv-vcf1 ref :interval-file target-bed)
        vcf-itree (parse-vcf-sv sv-vcf1 ref :out-format :itree)]
    (-> vcf-list first :start-ci) => 6066065
    (-> vcf-itree (get-itree-overlap "22" 15883520 15883620) first :end-ci) => 15883626
    (count vcf-by-region) => 1
    (with-open [vcf-source1 (get-vcf-source sv-vcf1 ref)
                vcf-source2 (get-vcf-source sv-vcf2 ref)]
      (doall (map #(get-sv-type % params) (parse-vcf vcf-source1))) =>
      (concat [:INS] (repeat 6 :BND)
              [nil :DEL :INS :DEL :DUP :INV :INS])
      (doall (map #(get-sv-type % params) (parse-vcf vcf-source2))) =>
      [:DUP :BND :BND :INS :CNV :DEL :INV])))

(facts "Compare structural variation calls from two inputs."
  (compare-sv "Test" {:name "sv1000g" :file sv-vcf1}
              {:name "svIll" :file sv-vcf2} ref) => (contains sv-out)
  (compare-sv "Test" {:name "sv1" :file sv-vcf1}
              {:name "sv2" :file sv-vcf1} ref) => (contains sv-out2))

(facts "Combine indels from different calling methodologies that overlap."
  (-> (compare-sv "Test" {:name "svindfb" :file indel-vcf1} {:name "svindgatk" :file indel-vcf2}
                  ref :params {:min-indel 2 :default-cis [[100 10]]})
      :sv-concordant
      (get-vcf-source ref)
      parse-vcf
      count) => 18)
