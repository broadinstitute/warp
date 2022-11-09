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
