(ns bcbio.variation.test.complex
  "Tests for dealing with more complex variations: structural
  variations and MNPs"
  (:use [midje.sweet]
        [bcbio.variation.combine :only [full-prep-vcf]]
        [bcbio.variation.compare :only [variant-comparison-from-config]]
        [bcbio.variation.complex]
        [bcbio.variation.normalize]
        [bcbio.variation.structural]
        [bcbio.variation.variantcontext])
  (:require [me.raynes.fs :as fs]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.variation.utils.svmerge :as svmerge]))

(background
 (around :facts
         (let [data-dir (str (fs/file "." "test" "data"))
               ref (str (fs/file data-dir "GRCh37.fa"))
               target-bed (str (fs/file data-dir "target-regions.bed"))
               sv-vcf1 (str (fs/file data-dir "sv-1000g.vcf"))
               sv-vcf2 (str (fs/file data-dir "sv-illumina.vcf"))
               multi-vcf (str (fs/file data-dir "1000genome-multi.vcf"))
               multi-out (fsp/add-file-part multi-vcf "nomnp")
               sv-out {:sv-concordant
                       (str (fs/file data-dir "sv-sv1000g-svIll-svconcordance.vcf"))
                       :sv-sv1000g-discordant
                       (str (fs/file data-dir "sv-sv1000g-svIll-svdiscordance.vcf"))
                       :sv-svIll-discordant
                       (str (fs/file data-dir "sv-svIll-sv1000g-svdiscordance.vcf"))}
               sv-out2
               {:sv-concordant (str (fs/file data-dir "sv-sv1-sv2-svconcordance.vcf"))
                :sv-sv1-discordant (str (fs/file data-dir "sv-sv1-sv2-svdiscordance.vcf"))
                :sv-sv2-discordant (str (fs/file data-dir "sv-sv2-sv1-svdiscordance.vcf"))}
               mnp-vcf (str (fs/file data-dir "freebayes-calls-indels.vcf"))
               cindel-vcf (str (fs/file data-dir "freebayes-calls-complexindels.vcf"))
               cindel-out (fsp/add-file-part cindel-vcf "nomnp")
               cindel-extras (map #(fsp/add-file-part cindel-vcf %)
                                  ["leftalign" "worknomnp" "worknomnp-leftalign"])
               indel-vcf1 (str (fs/file data-dir "sv-indels-fb.vcf"))
               indel-vcf2 (str (fs/file data-dir "sv-indels-gatk.vcf"))
               indel-out (str (fs/file data-dir "Test-svindfb-svindgatk-svconcordance.vcf"))
               nomnp-out (fsp/add-file-part mnp-vcf "nomnp")
               fullprep-out (fsp/add-file-part mnp-vcf "fullprep")
               headerfix-out (fsp/add-file-part mnp-vcf "samplefix")
               ;; SV testing
               presv-orig (str (fs/file data-dir "svs" "reference-calls.vcf"))
               sv-calls (str (fs/file data-dir "svs" "reference-calls-sv.vcf"))
               presv-regions (str (fs/file data-dir "svs" "reference-regions.bed"))
               sv-out-calls (fsp/add-file-part presv-orig "wsvs")
               sv-out-regions (fsp/add-file-part presv-regions "wsvs")
               sv-workdir (str (fs/file data-dir "svs" "work"))
               presv-extras [(fsp/add-file-part sv-calls "prep")
                             (fsp/add-file-part presv-orig "safesv")]
               params {:max-indel 100}]
           (doseq [x (concat [nomnp-out indel-out cindel-out headerfix-out fullprep-out
                              multi-out sv-workdir sv-out-calls sv-out-regions]
                             presv-extras cindel-extras (vals sv-out) (vals sv-out2))]
             (fsp/remove-path x))
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
    (with-open [vcf-iter1 (get-vcf-iterator sv-vcf1 ref)
                vcf-iter2 (get-vcf-iterator sv-vcf2 ref)]
      (doall (map #(get-sv-type % params) (parse-vcf vcf-iter1))) =>
      (concat [:INS] (repeat 6 :BND)
              [nil nil nil :DEL :INS :DEL :DUP :INV :INS])
      (doall (map #(get-sv-type % params) (parse-vcf vcf-iter2))) =>
      [:DUP :UNASSEMBLED_EVENT :BND :BND :INS nil nil :CNV :DEL :INV])))

(facts "Compare structural variation calls from two inputs."
  (compare-sv {:name "sv1000g" :file sv-vcf1}
              {:name "svIll" :file sv-vcf2} ref
              :params {:default-cis [[200 10]] :max-indel 30}) => (contains sv-out)
  (compare-sv {:name "sv1" :file sv-vcf1}
              {:name "sv2" :file sv-vcf1} ref) => (contains sv-out2))

(facts "Compare structural variants considering overlapping smaller variants"
  (let [cfile (str (fs/file data-dir "svs" "sv-grading.yaml"))]
    (svmerge/into-calls presv-orig presv-regions sv-calls ref) =>
      {:calls sv-out-calls :regions sv-out-regions}
    (let [sum (-> (variant-comparison-from-config cfile) first :summary)]
      (get-in sum [:sv :sv-reference-discordant :total]) => 1
      (get-in sum [:concordant :total]) => 2
      (get-in sum [:discordant2 :total]) => 0)))

(facts "Combine indels from different calling methodologies that overlap."
  (-> (compare-sv {:name "svindfb" :file indel-vcf1} {:name "svindgatk" :file indel-vcf2}
                  ref :params {:max-indel 2 :default-cis [[100 10]]})
      :sv-concordant
      (get-vcf-iterator ref)
      parse-vcf
      count) => 18)

(facts "Full preparation pipeline for variant files"
  (full-prep-vcf mnp-vcf ref) => fullprep-out)

(facts "Fix VCF header problems"
  (fix-vcf-sample mnp-vcf "Test1" ref) => headerfix-out)

(facts "Handle preparation of multi-sample input files"
  (normalize-variants multi-vcf ref) => multi-out)

(facts "Normalize problematic input VCF files, handling special cases"
  (let [inp-vcf (str (fs/file data-dir "phasing-comparison-needprep.vcf"))
        outp-vcf (fsp/add-file-part inp-vcf "fullprep")]
    (fsp/remove-path outp-vcf)
    (full-prep-vcf inp-vcf ref :keep-ref true) => outp-vcf))
