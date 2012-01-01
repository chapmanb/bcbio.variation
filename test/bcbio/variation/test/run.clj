;; Test code for idempotent file running

(ns bcbio.variation.test.run
  (:use [midje.sweet]
        [bcbio.run.itx])
  (:require [fs.core :as fs]))

(fact "Idempotent processing of files"
  (needs-run? "project.clj") => false
  (needs-run? "noexist.txt") => true
  (needs-run? "project.clj" "README.md") => false
  (needs-run? ["project.clj" "README.md"]) => false
  (needs-run? "project.clj" "noexist.txt") => true)

(facts "Manipulating file paths"
  (add-file-part "test.txt" "add") => "test-add.txt"
  (add-file-part "/full/test.txt" "new") => "/full/test-new.txt"
  (file-root "/full/test.txt") => "/full/test")
