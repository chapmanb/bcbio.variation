(ns bcbio.variation.test.multiple
  "Tests for post-comparison processing of multiple samples."
  (:use [midje.sweet]
        [bcbio.variation.compare :only [variant-comparison-from-config]]
        [bcbio.variation.config :only [load-config]]
        [bcbio.variation.filter.trusted]
        [bcbio.variation.filter.train :only [extract-train-cases]]
        [bcbio.variation.multiple]
        [bcbio.variation.variantcontext :exclude [-main]])
  (:require [fs.core :as fs]
            [bcbio.run.itx :as itx]))

(defn- prep-variant-comparison* [out-dir config-file]
  (itx/remove-path out-dir)
  (variant-comparison-from-config config-file))
(def prep-variant-comparison (memoize prep-variant-comparison*))

(background
 (around :facts
         (let [config-file (str (fs/file "." "config" "method-comparison.yaml"))
               config (load-config config-file)
               out-dir (str (fs/file (get-in config [:dir :prep]) "multiple"))
               union-file (str (fs/file (get-in config [:dir :prep]) "multiple"
                                        "Test1-multiall-fullcombine-gatk-annotated.vcf"))
               train-dir (str (fs/file (get-in config [:dir :prep]) "train"))
               trusted-out (itx/add-file-part union-file "trusted")
               cmps (prep-variant-comparison out-dir config-file)]
           (doseq [x [trusted-out train-dir]]
             (itx/remove-path x))
           ?form)))

(defn get-out-files [out-dir x ext]
  {:true-positives
   (str (fs/file out-dir (format "Test1-multiall-fullcombine-Intersection%s.vcf" ext)))
   :false-negatives
   (str (fs/file out-dir (format "Test1-multiall-no%s-fullcombine-%s%s.vcf" x x ext)))
   :false-positives
   (str (fs/file out-dir (format "Test1-dis%s-fullcombine-Intersection-shared.vcf" x)))
   :target-overlaps
   (str (fs/file out-dir (format "Test1-multiall-fullcombine-%s%s.vcf" x ext)))})
   
(facts "Handle multiple variant approach comparisons."
  (multiple-overlap-analysis cmps config "cg") => (get-out-files out-dir "cg" "")
  (multiple-overlap-analysis cmps config "gatk") => (get-out-files out-dir "gatk" "-annotated"))

(facts "Prepare trusted variant file"
  (get-trusted-variants cmps "gatk"
                        {:total 3 :technology 2}
                        (-> config :experiments first) config) => trusted-out)

(defn- get-train-out
  [dname base]
  (str (fs/file dname (format "Test1-freebayes-potential-%s.vcf" base))))

(facts "Identify true/false positives/negatives for training"
  (let [train-cfg {:target "freebayes"
                   :cmps ["gatk" "cg"]}
        exp (-> config :experiments first)]
    (extract-train-cases cmps train-cfg exp config)) => {:fps (get-train-out train-dir "fps")
                                                         :fns (get-train-out train-dir "fns")})

(facts "Prepare trusted variant file"
  (with-open [vcf-iter (get-vcf-iterator union-file (-> config :experiments first :ref))]
    (let [vc (first (parse-vcf vcf-iter))
          calls (-> config :experiments first :calls)
          data (variant-set-metadata vc calls)]
      (-> data :total count) => 3
      (-> data :technology) => #{"illumina" "cg"}
      (is-trusted-variant? vc {:total 4} calls) => nil
      (is-trusted-variant? vc {:total 3} calls) => true
      (is-trusted-variant? vc {:technology 3} calls) => nil
      (is-trusted-variant? vc {:technology 2} calls) => true
      (is-trusted-variant? vc {:technology 0.99} calls) => true
      (is-trusted-variant? vc {:total 4 :technology 3} calls) => nil 
      (is-trusted-variant? vc {:total 4 :technology 2} calls) => true)))