(ns bcbio.variation.web.shared
  "Shared functionality useful across multiple pages.")
  
(def ^{:doc "Web configuration, loaded from input YAML file"}
  web-config (atom nil))
