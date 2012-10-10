(ns bcbio.variation.remote.core
  "Top level API for working with remote filestores."
  (:require [bcbio.variation.remote.file :as file]
            [bcbio.variation.remote.client :as client]))

(def get-client client/get-client)
(def list-dirs file/list-dirs)
(def list-files file/list-files)
(def get-file file/get-file)
(def put-file file/put-file)