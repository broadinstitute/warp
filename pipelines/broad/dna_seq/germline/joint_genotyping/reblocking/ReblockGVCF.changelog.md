# 2.1.4
2022-11-04 (Date of Last Commit)

* Updated GATK verison to 4.3.0.0

# 2.1.3
2022-06-21 (Date of Last Commit)

* Changed QC.CheckFingerprint to QC.CheckFingerprintTask to avoid a naming conflict in the update scala tests, no effect on this pipeline

# 2.1.2
2022-06-01 (Date of Last Commit)

* Added inputs to the GenotypeGVCFs task to support the UltimaGenomicsJointGenotyping.wdl

# 2.1.1
2022-04-21 (Date of Last Commit)

* Fixed path to docker image in GermlineVariantDiscovery.wdl

# 2.1.0
2022-04-14 (Date of Last Commit)

* Updated to Picard version 2.26.10 and GATK version 4.2.6.1 to address log4j vulnerabilities
    * Reblocking fix to merge sites with missing DP into adjacent ref blocks
    
# 2.0.5
2022-03-24 (Date of Last Commit)

* Task wdls used by the ReblockGVCF pipeline were updated with changes that don't affect the ReblockGVCF pipeline itself

# 2.0.4
2022-02-01 (Date of Last Commit)

* Increase disk for Reblock Task

# 2.0.3
2022-01-14 (Date of Last Commit)

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