## 0.0.9 (in progress)

- Move to GATK 2.5-2 from GATK 2.3 lite. Required using new GATK framework and
  porting of used walkers (LeftAlignVariants) and annotators from GATK-lite code
  base. Removal of GATK 2.5 specific functionality which is no longer MIT licensed:
  recalling with UnifiedGenotyper and Variant Quality Score Recalibration.

- Improve normalization of complex variants (MNPs and indels) with handling of
  additional tricky cases. Thanks to David Mittelman and Jason Wang for example
  cases.

- Work towards slimming down main distribution to remove larger external dependencies.

## 0.0.8 (6 May 2013)

- Provide detailed breakdown of variant evaluation and grading to interface with
  bcbio-nextgen.

- Fixes for normalization to handle tricky structural variation cases. Thanks to
  David Jenkins for examples.

- Provide additional variant pre-cleaning steps to handle variant inputs not
  parseable by GATK.
