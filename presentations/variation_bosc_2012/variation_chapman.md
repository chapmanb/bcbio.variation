% Toolkit for variation comparison and analysis
% Brad Chapman, Bioinformatics Core at Harvard School of Public Health
% Bioinformatics Open Source Conference, 13 July 2012

# Variation

    1	1156131	rs2887286	T	C	1714.07	PASS	
    AB=0;ABP=0;AC=1;AC1=2;AF=1.0;AF1=1;AN=1;AO=56;
    BVAR;CIGAR=1X;DB;DP=59;DP4=0,0,27,26;DPRA=0;Dels=0.00;
    EPP=6.88793;EPPR=0;FQ=-187;FS=0.000;GC=55.45;HRun=0;
    HaplotypeScore=0.0000;LEN=1;MEANALT=2;MQ=70.00;MQ0=0;
    MQM=70;MQMR=0;NBQ=27.76;NS=1;NUMALT=1;ODDS=308.76;
    PAIRED=0.928571;QD=30.61;RO=0;RPP=12.937;RUN=1;
    SAP=3.16541;SNPEFF_EFFECT=INTRON;SNPEFF_GENE_NAME=SDF4;
    SNPEFF_FUNCTIONAL_CLASS=NONE;SNPEFF_IMPACT=MODIFIER;
    SNPEFF_GENE_BIOTYPE=protein_coding;
    SNPEFF_TRANSCRIPT_ID=ENST00000263741;
    Samples=NA19239_illumina;TYPE=snp;VDB=0.0416;
    VQSLOD=16.7453;XAI=0.00018797;XAM=0.00154116
    GT:AO:DP:GL:GQ:QA:QR:RO
    1:56:58:-141.12,-3.61:99:1527:0:0
    
---

\includegraphics[width=1.0\textwidth]{images/variation_abstraction}

# Answer biological questions; help people

\includegraphics[width=1.0\textwidth]{images/cmt_goals}

\url{http://cmtproject.blogspot.com/}
    
# Solutions

Comparisons

- Multiple callers: GATK, FreeBayes, samtools
- Multiple technologies: Illumina, SOLiD, IonTorrent

Identify real variants

- Summarize associated metrics
- Remove false positives

Scale

- Millions of variants
- Thousands of samples

---

\includegraphics[width=0.7\textwidth]{images/xprize_logo.png}
\vspace{1em}
\includegraphics[width=0.6\textwidth]{images/edgebio_logo.png}

\url{http://genomics.xprize.org/}

# Clinical grade genome

- 98 percent genome coverage
- 1 error per million bases (SNPs + small indels)
- Full haplotype phasing
- Structural variations

# Sequencing for patients
\includegraphics[width=1.0\textwidth]{images/bertrand_intro}
\vspace{0.1em}
\url{http://matt.might.net/articles/my-sons-killer/}

# Technology overview

- Clojure
  \vspace{0.1em}
  \url{http://clojure.org/}
- Genome Analysis Toolkit
  \vspace{0.1em}
  \url{http://www.broadinstitute.org/gsa/wiki/index.php/Home_Page}
- GenomeSpace
  \vspace{0.1em}
  \url{http://www.genomespace.org/}  

\url{https://github.com/chapmanb/bcbio.variation}

# Why Clojure?

- Hosted on Java Virtual Machine
    - Interoperability with existing libraries: GATK, GenomeSpace
    - Wonderful build, deployment and testing tools
- Functional and immutable
    - Easier to write correct code
    - Small functions: concise and refactorable
- Community
    - Smart people working on hard problems
- Ecosystem
    - Multiple backends: ClojureScript

# Example analysis pipeline

- Variant files from two different callers
    - GATK UnifiedGenotyper
    - FreeBayes: \url{https://github.com/ekg/freebayes}
- Compare, identifying:
    - Identical variants
    - Different variants in each caller
    - Metrics to help discriminate

# High level description

    experiments:
      - sample: Test1
        ref: test/data/GRCh37.fa
        intervals: test/data/target-regions.bed
        align: test/data/aligned-reads.bam
        calls:
          - name: gatk
            file: test/data/gatk-calls.vcf
          - name: freebayes
            file: test/data/freebayes-calls.vcf
            prep: true
            annotate: true
            filters:
              - QD < 2.0
              - MQRankSum < -12.5

# Simple to run

With a new blank machine, get Java and automated Clojure build tool:

    $ sudo apt-get install openjdk-7-jdk openjdk-7-jre
    $ wget https://raw.github.com/technomancy/leiningen/\
           preview/bin/lein
    $ chmod 755 lein && sudo mv lein /usr/local/bin
    
Get code and dependencies, then run:
    
    $ git pull git://github.com/chapmanb/bcbio.variation.git
    $ cd bcbio.variation
    $ lein deps
    $ lein variant-compare comparison_description.yaml
   
# What happened?

\includegraphics[width=1.0\textwidth]{images/comparison_workflow}

# Establishing true variants

X Prize: haploid gold standard genome

- Public genomes from HapMap/1000 genomes
    - NA12878: Caucasian female from Utah.
    - NA19239: Yoruba male from Ibadan, Nigeria.
    
- Multiple technologies
    - Illumina
    - SOLiD
    - IonTorrent
    
- Multiple callers
    - GATK
    - FreeBayes
    - samtools mpileup

# True variant workflow -- comparisons

\includegraphics[width=1.0\textwidth]{images/prep_reference_input}

# True variant workflow -- finalize

\includegraphics[width=1.0\textwidth]{images/prep_reference_recal}

# Comparison web interface

\includegraphics[width=1.0\textwidth]{images/web_compare}

\url{https://validationprotocol.org/}

# Summary web interface

\includegraphics[width=1.0\textwidth]{images/web_summary}

# Next web steps -- metrics

\includegraphics[width=1.0\textwidth]{images/web_metrics}

\url{http://keminglabs.com/}

---

\includegraphics[width=1.0\textwidth]{images/variation_abstraction}

# Answering biological questions

- Establish set of true variants
    - Understand boundaries of certainty
    - Make patient decisions
    
- Comparison architecture
    - Cancer: tumor/normal pairs
    - Mendelian inherited diseases: father/mother/child
    
- Annotate with known data
    - Gemini: \url{https://github.com/arq5x/gemini}
    - 1000 genomes frequencies
    - Mappability
    - Clinical information
