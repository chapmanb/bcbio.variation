(ns bcbio.run.parallel
  "Utilities for running code in parallel."
  (:require [clojure.core.reducers :as r]))

(defn rmap
  "Reducer based parallel map, with flexible core usage and chunking.
   http://www.thebusby.com/2012/07/tips-tricks-with-clojure-reducers.html
   https://groups.google.com/d/msg/clojure/asaLNnM9v74/-t-2ZlCN5P4J
   https://groups.google.com/d/msg/clojure/oWyDP1JGzwc/5oeYqEHHOTAJ"
  ([f coll cores chunk-size]
     (if (= 1 cores)
       (map f coll)
       (do
         (alter-var-root #'r/pool (constantly (future (java.util.concurrent.ForkJoinPool. (int cores)))))
         (r/fold chunk-size r/cat r/append! (r/map f (vec coll))))))
  ([f coll cores]
     (rmap f coll cores 1))
  ([f coll]
     (rmap f coll 1 1)))
