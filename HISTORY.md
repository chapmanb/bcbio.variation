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
