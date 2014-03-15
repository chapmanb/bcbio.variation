## 0.1.5 (15 March 2014)

- Move to MIT licensed GATK 3.0 framework.
- Support bgzipped inputs for variant assessment. Thanks to Severine Catreux.
- Support lightweight loading options for gemini integration to avoid large load
  times with new gene tables in gemini 0.6.5.
- Avoid running GATK VariantEval which causes intermittent java core dumps.
- Avoid re-runs of callable regions when already prepared using non-GATK chanjo-based
  methods.

## 0.1.4 (5 March 2014)

- Allow bgzipped/tabix inputs to validation.
- Improve representation of ensemble variants ensuring standard expected keys
  are present when available in an individual genotype call. Thanks to Shalabh
  Suman.

## 0.1.3 (15 January 2014)

- Update dependencies to the GATK 2.8.1 MIT licensed framework and associated
  Picard, Tribble and variant libraries. Silence phone home events from library.
- Avoid errors on converting hg19 to GRCh37 where hg19 variants contain hg19 hap
  contigs with no equivalent in GRCh37. It now drops these variants instead
  of generating an error. Thanks to Severine Catreux.
- Ensure BED files passed to evaluations are converted to reference material
  coordinates when inputs differ (hg19->GRCh37). Thanks to Severine Catreux.
- Avoid issues with running LeftAlignVariants on indels with END tags. Thanks to
  Justin Johnson.
- Move to external bcbio.run tool to help abstract out some core functionality
  useful in other contexts.
- Additional flexibility for Illumina to valid VCF preparation. Allow ignoring
  SVs or other file types.

## 0.1.2 (4 December 2013)

- Handle metrics evaluation with GATK where input calls have haploid chromosomes,
  in the case of mitochondrial and chrY for males.
- Check configuration input keys for comparisons to try and provide useful error
  messages for unexpected keys.
- Provide `variant-utils sort-vcf`, which sorts a VCF file to match a reference
  genome file.

## 0.1.1 (20 October 2013)

- Better defaults for selecting false positives during Ensemble calling to
  handle common 3-caller use cases.
- Improve output for Ensemble based calling, providing useful `set` INFO
  attribute. Thanks to Shalabh Suman for discussions.
- Increase length requirement for exact comparisons to 100bp to account for
  improved resolution at longer read sizes.
- Correctly handle overlapping structural variations of same type.
- Improve sorting of BED files to be faster and less memory dependent. Thanks to
  Zhengqiu Cai for problem report.
- Move to GATK 2.7 framework. Eliminated dependencies to reduce size of final
  jar.
- Custom filtering with complex attribute comparisons, allowing usage as a
  commandline program.
- Move threshold for low-coverage in discordant identification to 13 based on
  http://www.ncbi.nlm.nih.gov/pubmed/23773188

## 0.1.0 (25 August 2013)

- Provide `ensemble` run target to allow easier normalization of inputs from
  multiple variant files.
- Generalize filtering and combining metrics for multi-sample inputs. Allows
  ensemble calls on jointly called VCFs.
- Speed improvements for normalizing input VCFs in preparation for comparison or
  downstream manipulation.
- Handle merging of Illumina structural variant calls (SVs.vcf) during
  preparation of Illumina VCFs to GATK friendly VCFs.
- Bug fixes for removing variants that don't match reference as part of the
  variant preparation. Now handle hg19 -> GRCh37 chromosome swaps and
  multibase deletions.
- Bug fix for Neighboring Base Quality annotation metric for regions with no
  quality scores. Thanks to Zhengqiu Cai.

## 0.0.9 (22 June 2013)

- Move to GATK 2.5-2 from GATK 2.3 lite. Required using new GATK framework and
  porting of used walkers (LeftAlignVariants) and annotators from GATK-lite code
  base. Removal of GATK 2.5 specific functionality which is no longer MIT licensed:
  recalling with UnifiedGenotyper and Variant Quality Score Recalibration.

- Improve normalization of complex variants (MNPs and indels) with handling of
  additional tricky cases. Thanks to David Mittelman and Jason Wang for example
  cases.

- Work towards slimming down main distribution to remove larger external
  dependencies. Next releases will move external functionality into own packages.

## 0.0.8 (6 May 2013)

- Provide detailed breakdown of variant evaluation and grading to interface with
  bcbio-nextgen.

- Fixes for normalization to handle tricky structural variation cases. Thanks to
  David Jenkins for examples.

- Provide additional variant pre-cleaning steps to handle variant inputs not
  parseable by GATK.
