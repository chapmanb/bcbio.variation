# bcbio.variation

Use the [Genome Analysis Toolkit (GATK)][1] to analyze variant data.
This is a Clojure API to parse and analyze [VCF files][2]. It supports scoring
for the [Archon Genomic X PRIZE competition][5] but is also a general framework
for variant file comparison.

For more details, see the [code documentation][3] and a [presentation overview of the project][4].

[1]: http://www.broadinstitute.org/gsa/wiki/index.php/The_Genome_Analysis_Toolkit
[2]: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40
[3]: http://chapmanb.github.com/bcbio.variation
[4]: http://chapmanb.github.com/bcbio.variation/presentations/gatk_clojure.pdf
[5]: http://genomics.xprize.org/

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

[u1]: https://github.com/technomancy/leiningen
[u2]: http://en.wikipedia.org/wiki/YAML

## License

The code is freely available under the [MIT license][l1].

[l1]: http://www.opensource.org/licenses/mit-license.html
