# bcbio.variation

A Clojure interface to the [Genome Analysis Toolkit (GATK)][1] to analyze
variant data in [VCF files][2]. It supports scoring for the
[Archon Genomic X PRIZE competition][5] but is also a general framework
for variant file comparison.

* [presentation overview of the project][4]
* [howto description of interfacing with GATK][6]
* [code documentation][3]

[![Build Status](https://secure.travis-ci.org/chapmanb/bcbio.variation.png)](http://travis-ci.org/chapmanb/bcbio.variation)

[1]: http://www.broadinstitute.org/gsa/wiki/index.php/The_Genome_Analysis_Toolkit
[2]: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40
[3]: http://chapmanb.github.com/bcbio.variation
[4]: http://chapmanb.github.com/bcbio.variation/presentations/gatk_clojure.pdf
[5]: http://genomics.xprize.org/
[6]: http://bcbio.wordpress.com/2012/03/04/extending-the-gatk-for-custom-variant-comparisons-using-clojure/

## Usage

### Setup

Requires Java 1.6 and [Leiningen][u1].

    $ lein deps

### Generate summary of concordance between variant calls

A [YAML configuration file][u2] specifies the variant files for
comparison. The project contains example configuration and associated
variant files that demonstrate the features of the library.

An example of scoring a phased diploid genome against a haploid reference
genome:
    
    $ lein run :compare config/reference-grading.yaml


An example of assessing variant calls produced by different calling algorithms:

    $ lein run :compare config/method-comparison.yaml

### Web interface

A web interface automates the process of preparing configuration files and
running a variant comparison:
    
    $ lein run :web config/web-processing.yaml

### Run GATK walker for variant statistics

    $ lein uberjar
    $ java -jar bcbio.variation-0.0.1-SNAPSHOT-standalone.jar -T VcfSimpleStats
      -R test/data/GRCh37.fa --variant test/data/gatk-calls.vcf --out test.png

### Run custom GATK annotator

    $ lein uberjar
    $ java -jar bcbio.variation-0.0.1-SNAPSHOT-standalone.jar -T VariantAnnotator
       -A MeanNeighboringBaseQuality -R test/data/GRCh37.fa -I test/data/aligned-reads.bam
       --variant test/data/gatk-calls.vcf -o annotated-file.vcf

## Configuration file

A YAML configuration file defines targets for comparison processing. Two example
files for [reference grading][u4] and [comparison of calling methods][u3]
provide example starting points and details on available options are below:

    dir:
      base: Base directory to allow use of relative paths (optional).
      out: Working directory to write output.
      prep: Prep directory where files will be pre-processed.
    experiments: # one or more experiments
     - sample: Name of current sample.
       ref: Reference genome in FASTA format.
       intervals: Intervals to process in BED format (optional).
       align: Alignments for all calls in BAM format (optional).
       summary-level: Amount of summary information to provide,
                      [full,quick] (default:full)
       approach: Type of comparison to do [compare,grade]. Default compare.
       calls: # two or more calls to compare
         - name: Name of call type
           file: One or more input files in VCF format
           align: Alignment for specific call in BAM format (optional).
           ref: Reference genome if different than experiment ref (optional)
           intervals: Genome intervals to process in BED format (optional).
           refcalls: Add reference calls if has alignment info (boolean; default false).
           annotate: Annotate calls with GATK annotations (boolean; default false).
           normalize: Normalize MNPs and indels (boolean: default true).
           prep: Prep with in-order chromosomes and sample names (boolean; default false).
           preclean: Remove problematic characters from input VCFs (boolean; default false). 
           remove-refcalls: Remove reference, non-variant calls. (boolean; default false). 

[u1]: https://github.com/technomancy/leiningen
[u2]: http://en.wikipedia.org/wiki/YAML
[u3]: https://github.com/chapmanb/bcbio.variation/blob/master/config/method-comparison.yaml
[u4]: https://github.com/chapmanb/bcbio.variation/blob/master/config/reference-grading.yaml

## License

The code is freely available under the [MIT license][l1].

[l1]: http://www.opensource.org/licenses/mit-license.html
