# 2.4.4
2026-01-21 (Date of Last Commit)

* Moved inputs into new Google buckets. This change does not affect the outputs of the pipeline

# 2.4.3
2025-10-09 (Date of Last Commit)

* Modified the ReblockGVCF.wdl to use bash's basename instead of WDL's basename; this does not affect the outputs of this pipeline

# 2.4.2
2025-08-11 (Date of Last Commit)

* Reorganized pipelines into the wdl directory

# 2.4.1
2025-02-21 (Date of Last Commit)

* Updated HaplotypeCaller_GATK4_VCF to use MEM_SIZE and MEM_UNIT; this does not affect the outputs of this pipeline

# 2.4.0
2024-12-05 (Date of Last Commit)

* Updated output names for ReblockGVCF workflow from output_vcf and output_vcf_index to reblocked_gvcf and reblocked_gvcf_index respectively

# 2.3.2
2024-11-04 (Date of Last Commit)

* Updated to GATK version 4.6.1.0

# 2.3.1
2024-10-28 (Date of Last Commit)

* Updated GATK for Validate Variants, which reduces the memory requirements for the task when an interval list is not provided

# 2.3.0
2024-09-06 (Date of Last Commit)

* Updated to GATK version 4.6.0.0. Some expected minor differences around low quality sites (GQ0 genotypes or no-calls)

# 2.2.1
2024-06-12 (Date of Last Commit)

* ValidateVcfs is more robust to larger inputs

# 2.2.0
2024-07-09 (Date of Last Commit)

* Updated ReblockGVCF.wdl to run in Azure. cloud_provider is a new, required input 

# 2.1.13
2024-07-01 (Date of Last Commit)

* CalculateReadGroupChecksum requires more memory and disk; this does not affect this pipeline

# 2.1.12
2024-03-26 (Date of Last Commit)

* ValidateVcfs requires less memory when run without interval list. This is useful for WGS samples that were previously running out of memory

# 2.1.11
2023-12-18 (Date of Last Commit)

* Updated to GATK version 4.5.0.0. Header documentation change for RAW_GT_COUNT annotation

# 2.1.10
2023-12-14 (Date of Last Commit)

* Updated GATK for Reblock task to version 4.5.0.0
* Added options to Reblock task to remove annotations and move filters to genotype level

# 2.1.9
2023-12-08 (Date of Last Commit)

* ValidateVcfs now has optional memory parameter

# 2.1.8
2023-11-29 (Date of Last Commit)

* ValidateVcfs can now take in VCF as calling_interval_list that is in a separate location from its index

# 2.1.7
2023-11-21 (Date of Last Commit)

* Fixes bug so now ReblockGVCFs can take in GVCFs that are not in the same location as their index file

# 2.1.6
2023-09-18 (Date of Last Commit)

* ReblockGVCFs can now take in GVCFs that are not in the same location as their index file

# 2.1.5
2023-03-20 (Date of Last Commit)

* CheckFingerprint can allow LOD 0

# 2.1.4
2023-11-04 (Date of Last Commit)

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