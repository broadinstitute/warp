# 2.2.7
2025-10-09 (Date of Last Commit)

* Modified the ReblockGVCF.wdl to use bash's basename instead of WDL's basename; this does not affect the outputs of this pipeline

# 2.2.6
2025-08-11 (Date of Last Commit)

* Reorganized pipelines into the wdl directory

# 2.2.5
2025-02-21 (Date of Last Commit)

* Updated HaplotypeCaller_GATK4_VCF to use MEM_SIZE and MEM_UNIT; this does not affect the outputs of this pipeline

# 2.2.4
2024-11-04 (Date of Last Commit)

* Updated to GATK version 4.6.1.0

# 2.2.3
2024-10-28 (Date of Last Commit)

* Updated the docker in the ValidateVCF task; this does not affect this pipeline

# 2.2.2
2024-09-06 (Date of Last Commit)

* Updated to GATK version 4.6.0.0

# 2.2.1
2024-06-12 (Date of Last Commit)

* ValidateVcfs is more robust to larger inputs

# 2.2.0
2024-07-09 (Date of Last Commit)

* Updated tasks GermlineVariantDiscovery.wdl and QC.wdl to allow multi-cloud dockers. cloud_provider is a new, required input.

# 2.1.19
2024-07-01 (Date of Last Commit)

* CalculateReadGroupChecksum requires more memory and disk; this does not affect this pipeline

# 2.1.18
2024-03-26 (Date of Last Commit)

* ValidateVcfs requires less memory when run without interval list; this does not affect this pipeline

# 2.1.17
2023-12-18 (Date of Last Commit)

* Updated to GATK version 4.5.0.0. Header documentation change for RAW_GT_COUNT annotation.

# 2.1.16
2023-12-14 (Date of Last Commit)

* Updated GATK for Reblock task to version 4.5.0.0
* Added options to Reblock task to remove annotations and move filters to genotype level

# 2.1.15
2023-12-08 (Date of Last Commit)

* ValidateVcfs now has optional memory parameter; this does not affect this pipeline

# 2.1.14
2023-11-29 (Date of Last Commit)

* ValidateVcfs can now take in VCF as calling_interval_list that is in a separate location from its index; this does not affect this pipeline

# 2.1.13
2023-11-29 (Date of Last Commit)

* Fixed bug in ReblockGVCFs; this does not affect this pipeline.
* Reverted the VerifyBamID docker image back to the 2.1.10 VariantCalling pipeline version

# 2.1.12
2023-09-18 (Date of Last Commit)

* ReblockGVCFs can now take in GVCFs that are not in the same location as their index file, this update has no effect on this pipeline.

# 2.1.11
2023-08-23 (Date of Last Commit)

* Updated VerifyBamID docker image in BamProcessing.wdl to fix security vulnerabilities, this update has no effect on this pipeline.
* Updated the VCF validation step to only use \"--no-overlaps\" argument for reblocked vcfs

# 2.1.10
2023-03-20 (Date of Last Commit)

* CheckFingerprint can allow LOD 0

# 2.1.9
2022-11-04 (Date of Last Commit)

* Updated GATK verison to 4.3.0.0

# 2.1.8
2022-09-27 (Date of Last Commit)

* Removed task MakeOptionalOutputBam in Utilities.wdl, this update has no effect on this pipeline

# 2.1.7
2022-09-23 (Date of Last Commit)

* Updated Picard-Python Docker image in Utilities.wdl to fix vulnerabilities.

# 2.1.6
2022-07-15 (Date of Last Commit)

* Updated task MakeOptionalOutputBam in Utilities.wdl, this update has no effect on this pipeline

# 2.1.5
2022-07-12 (Date of Last Commit)

* Added additional_disk input to SortSam task in BamProcessing.wdl

# 2.1.4
2022-07-11 (Date of Last Commit)

* Added memory_multiplier and additional_disk inputs to GatherSortedBamFiles task in BamProcessing.wdl

# 2.1.3
2022-06-21 (Date of Last Commit)

* Changed QC.CheckFingerprint to QC.CheckFingerprintTask to avoid a naming conflict in the update scala tests, no effect on this pipeline

# 2.1.2
2022-06-01 (Date of Last Commit)

* Added a MakeOptionalOutputBam task to the Utilities.wdl tp support UltimaGenomicsWholeGenomeCramOnly.wdl
* Added inputs to the GenotypeGVCFs task to support the UltimaGenomicsJointGenotyping.wdl


# 2.1.1
2022-04-21 (Date of Last Commit)

* Fixed path to docker image in GermlineVariantDiscovery.wdl

# 2.1.0
2022-04-19 (Date of Last Commit)

* Updated to Picard version 2.26.10 and GATK version 4.2.6.1 to address log4j vulnerabilities
    * RAW_RankSum NaN to empty for NON_REF data 

# 2.0.6
2022-04-15 (Date of Last Commit)

* Updated task SortSam in BamProcessing.wdl to take an optional memory_multiplier

# 2.0.5
2022-04-04 (Date of Last Commit)

* Update task CopyFilesFromCloudToCloud in Utilities.wdl, this update has no effect on this pipeline

# 2.0.4
2022-03-24 (Date of Last Commit)

* Task wdls used by the VariantCalling pipeline were updated with changes that don't affect the VariantCalling pipeline itself

# 2.0.3
2022-02-01 (Date of Last Commit)

* Increased the disk space in Reblock task
* Increased the disk space in CalibrateDragstrModel task
* Addressed memory usage in CheckFingerprint task to allow sufficient headroom for the VM

# 2.0.2
2022-01-14 (Date of Last Commit)

* Task wdls used by VariantCalling were updated with changes that don't affect VariantCalling wdl

# 2.0.1
2021-12-9

* Updated broken dependency in VerifyBamID docker image

# 2.0.0
2021-11-15

* Added an optional step to reblock gVCFs, this step is included by default
    * The WholeGenomeGermlineSingleSample and ExomeGermlineSingleSample pipelines now output reblocked gVCFs by default. To skip reblocking, add '\"WholeGenomeGermlineSingleSample.BamToGvcf.skip_reblocking\": true' or '\"ExomeGermlineSingleSample.BamToGvcf.skip_reblocking\": true' to the appropriate inputs
* Added WGS plumbing tests for dragen_maximum_quality_mode and dragen_functional_equivalence_mode
* Added Xmx flag (maximum heap size) to all tasks with java commands
* Task wdls used by the VariantCalling pipeline were updated with changes that don't affect VariantCalling wdl

# 1.2.0
2021-10-18

* Added new optional workflow inputs to support the DRAGEN-GATK mode of the Whole Genome Germline Single Sample (WGS) workflow (read more in the [WGS Overview](https://broadinstitute.github.io/warp/docs/Pipelines/Whole_Genome_Germline_Single_Sample_Pipeline/README)). These include:
    * `Boolean run_dragen_mode_variant_calling = false`
    * `Boolean use_spanning_event_genotyping = true`
    * `Boolean use_dragen_hard_filtering = false`
    * `File? ref_str`
* Added a new task (DragenTasks) to support variant calling in DRAGEN mode
* Updated GATK to v4.2.2.0 for variant calling

# 1.1.0
2021-10-06

* Updated VerifyBamID to use AppSec base image
* Changed GoTC image to Samtools specific image in CramToUnmappedBams and Utilities
* Changed GoTC image to GATK specific image in GermlineVariantDiscovery

# 1.0.2
2021-09-22

* Updated Utilities.wdl task definitions to include a new ErrorWithMessage task that is NOT used in the VariantCalling pipeline.

# 1.0.1
2021-06-22

* Removed duplicate MarkDuplicatesSpark task from BamProcessing
* Removed duplicate Docker image from CheckPreValidation task in QC

# 1.0.0
2021-03-17

* Promoted VariantCalling to be a top-level workflow.
