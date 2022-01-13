# 2.0.3
2022-01-12

* Task wdls used by the ReblockGVCF pipeline were updated with changes that don't affect ReblockGVCF wdl

# 2.0.2
2021-11-15

* Task wdls used by the ReblockGVCF pipeline were updated with changes that don't affect ReblockGVCF wdl

# 2.0.1
2021-11-10

* Added a validation task to check output reblocked GVCFs for reference block overlaps
* Moved ReblockGVCF task to GermlineVariantDiscovery tasks wdl
* Added WGS plumbing tests for dragen_maximum_quality_mode and dragen_functional_equivalence_mode
* Added Xmx flag (maximum heap size) to all tasks with java commands

# 2.0.0
2021-08-17

Updated to ReblockGVCF in [GATK 4.2.2.0](https://github.com/broadinstitute/gatk/releases/tag/4.2.2.0).  Now output GVCFs: 
  *  Cover every position
  *  Do not contain overlapping reference blocks
  *  Have correct reference allele following trimmed deletions
Tool, task, and workflow now require a reference.

# 1.1.0
2020-09-25

* Updated GATK docker image for all tasks from [GATK 4.1.1.0](https://github.com/broadinstitute/gatk/releases/tag/4.1.1.0) to [GATK 4.1.8.0](https://github.com/broadinstitute/gatk/releases/tag/4.1.8.0).
    * Numerous bug fixes and improvements, including better VQSR reliability for allele-specific filtering.
    * Addition of new Scattered CrosscheckFingerprints to pipeline

# 1.0.1
2020-08-28

* Added 'allowNestedInputs: true' metadata parameter to wdl to support Cromwell version 48 and above

# 1.0
Initial release of the ReblockGVCF pipeline. This is a WDL 1.0 version of the [workflow run in Terra](https://portal.firecloud.org/?return=terra#methods/methodsDev/ReblockGVCF-gatk4_exomes_goodCompression/4) for joint-genotyping projects. 