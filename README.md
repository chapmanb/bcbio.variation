# bcbio.variation

Use the [Genome Analysis Toolkit (GATK)][1] to analyze variant data.
This is a Clojure API to parse and analyze [VCF files][2].
See the [documentation][3] for more details.

[1]: http://www.broadinstitute.org/gsa/wiki/index.php/The_Genome_Analysis_Toolkit
[2]: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40
[3]: http://chapmanb.github.com/bcbio.variation

## Usage

### Load required libraries

Requires Java 1.6 and [Leiningen][3].

    $ lein deps

### Generate summary of concordance between variant calls

    $ lein run :compare config/method-comparison.yaml

### Run GATK walker for variant statistics

    $ lein uberjar
    $ java -jar bcbio.variation-0.0.1-SNAPSHOT-standalone.jar -T VcfSimpleStats
      -R test/data/hg19.fa --variant test/data/gatk-calls.vcf --out test.png

[3]: https://github.com/technomancy/leiningen

## License

The code is freely available under the [MIT license][l1].

[l1]: http://www.opensource.org/licenses/mit-license.html
