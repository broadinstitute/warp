# 3.3.7
2026-01-21 (Date of Last Commit)

* Moved inputs into new Google buckets. This change does not affect the outputs of the pipeline

# 3.3.6
2025-10-09 (Date of Last Commit)

* Modified the ReblockGVCF.wdl to use bash's basename instead of WDL's basename; this does not affect the outputs of this pipeline

# 3.3.5
2025-08-11 (Date of Last Commit)

* Reorganized pipelines into the wdl directory

# 3.3.4
2025-02-21 (Date of Last Commit)

* Updated HaplotypeCaller_GATK4_VCF to use MEM_SIZE and MEM_UNIT; this does not affect the outputs of this pipeline

# 3.3.3
2024-11-04 (Date of Last Commit)

* Updated to GATK version 4.6.1.0

# 3.3.2
2024-10-28 (Date of Last Commit)

* Updated the docker in the ValidateVCF task; this does not affect this pipeline

# 3.3.1
2024-09-17 (Date of Last Commit)

* Updated DRAGEN aligner parameters to fix non-determinism

# 3.3.0
2024-09-06 (Date of Last Commit)

* Updated to GATK version 4.6.0.0
* Updated version drops some low quality sites from VCFs; if reblocking is enabled, the DP annotation in some ref blocks will change due to the change in HaplotypeCaller

# 3.2.2
2024-08-02 (Date of Last Commit)

* The ubuntu_16_0_4 docker image version was pinned instead of using the latest tag; this does not affect the outputs of the pipeline

# 3.2.1
2024-06-12 (Date of Last Commit)

* ValidateVcfs is more robust to larger inputs; this does not affect this pipeline

# 3.2.0
2024-07-09 (Date of Last Commit)

* Updated tasks GermlineVariantDiscovery.wdl and QC.wdl to allow multi-cloud dockers. cloud_provider is a new, required input.

# 3.1.21
2024-07-01 (Date of Last Commit)

* CalculateReadGroupChecksum requires more memory and disk

# 3.1.20
2024-03-26 (Date of Last Commit)

* ValidateVcfs requires less memory when run without interval list; this does not affect this pipeline

# 3.1.19
2023-12-18 (Date of Last Commit)

* Updated to GATK version 4.5.0.0.

# 3.1.18
2023-12-14 (Date of Last Commit)

* Updated GATK for Reblock task to version 4.5.0.0
* Added options to Reblock task to remove annotations and move filters to genotype level

# 3.1.17
2023-12-08 (Date of Last Commit)

* ValidateVcfs now has optional memory parameter; this does not affect this pipeline

# 3.1.16
2023-11-29 (Date of Last Commit)

* ValidateVcfs can now take in VCF as calling_interval_list that is in a separate location from its index; this does not affect this pipeline

# 3.1.15
2023-11-29 (Date of Last Commit)

* Fixed bug in ReblockGVCFs; this does not affect this pipeline.
* Reverted the VerifyBamID docker image back to the 3.1.11 WholeGenomeReprocessing pipeline version

# 3.1.14
2023-10-10 (Date of Last Commit)

* Removed the SumFloats task from SplitLargeReadGroup.wdl; this does not affect the outputs

# 3.1.13
2023-09-18 (Date of Last Commit)

* ReblockGVCFs can now take in GVCFs that are not in the same location as their index file, this update has no effect on this pipeline.

# 3.1.12
2023-08-23 (Date of Last Commit)

* Updated VerifyBamID docker image in BamProcessing.wdl to fix security vulnerabilities, this update has no effect on this pipeline.
* Updated the VCF validation step to only use \"--no-overlaps\" argument for reblocked vcfs

# 3.1.11
2023-03-20 (Date of Last Commit)

* CheckFingerprint can allow LOD 0

# 3.1.10
2022-11-04 (Date of Last Commit)

* Updated GATK verison to 4.3.0.0

# 3.1.9
2022-09-27 (Date of Last Commit)

* Removed task MakeOptionalOutputBam in Utilities.wdl, this update has no effect on this pipeline

# 3.1.8
2022-09-23 (Date of Last Commit)

* Updated Picard-Python Docker image in Utilities.wdl to fix vulnerabilities.

# 3.1.7
2022-07-15 (Date of Last Commit)

* Updated task MakeOptionalOutputBam in Utilities.wdl, this update has no effect on this pipeline

# 3.1.6
2022-07-12 (Date of Last Commit)

* Added additional_disk input to SortSam task in BamProcessing.wdl

# 3.1.5
2022-07-11 (Date of Last Commit)

* Added memory_multiplier and additional_disk inputs to GatherSortedBamFiles task in BamProcessing.wdl

# 3.1.4
2022-06-21 (Date of Last Commit)

* Changed QC.CheckFingerprint to QC.CheckFingerprintTask to avoid a naming conflict in the update scala tests, no effect on this pipeline

# 3.1.3
2022-06-01 (Date of Last Commit)

* Updated tasks in the QC.wdl and VariantCalling.wdl, this update has no effect on this pipeline 

# 3.1.2
2022-04-22 (Date of Last Commit)

* Updated task CopyFilesFromCloudToCloud in Utilities.wdl, this update has no effect on this pipeline

# 3.1.1
2022-04-21 (Date of Last Commit)

* Fixed path to docker image in GermlineVariantDiscovery.wdl

# 3.1.0
2022-04-19 (Date of Last Commit)

* Updated to Picard version 2.26.10 and GATK version 4.2.6.1 to address log4j vulnerabilities
    * The following metrics were added to alignment summary and readgroup alignment summary metrics:
        * AVG_POS_3PRIME_SOFTCLIP_LENGTH
        * MAD_READ_LENGTH
        * MAX_READ_LENGTH
        * MIN_READ_LENGTH
        * SD_READ_LENGTHMEDIAN_READ_LENGTH
    * Small differences observed in PCT_SOFTCLIP in alignment summary metrics due to a bug fix in the way PCT_SOFTCLIP is calculated
    * RAW_RankSum NaN to empty for NON_REF data 
    * Reblocking fix to merge sites with missing DP into adjacent ref blocks

# 3.0.7
2022-04-15 (Date of Last Commit)

* Updated task SortSam in BamProcessing.wdl to take an optional memory_multiplier

# 3.0.6
2022-04-04 (Date of Last Commit)

* Update task CopyFilesFromCloudToCloud in Utilities.wdl, this update has no effect on this pipeline

# 3.0.5
2022-03-24 (Date of Last Commit)

* The pipeline was modified to allow the read_length parameter to be overridden in the QC tasks CollectWGSMetrics and CollectRawWGSMetrics

# 3.0.4
2022-02-02 (Date of Last Commit)

* Changed dragmap base image from Centos to RockyLinux to comply with trivy scans

# 3.0.3
2022-02-01 (Date of Last Commit)

* Increased the disk space in Reblock task
* Increased the disk space in CalibrateDragstrModel task
* Addressed memory usage in CheckFingerprint task to allow sufficient headroom for the VM

# 3.0.2
2022-01-14 (Date of Last Commit)

* Increased the disk space in CalibrateDragstrModel task

# 3.0.1
2021-12-09
* Updated the base image for the Dragmap docker image
* Updated broken dependency in VerifyBamID docker image

# 3.0.0
2021-11-15

* Added an optional step to reblock gVCFs, this step is included by default
    * The WholeGenomeReprocessing pipeline now outputs reblocked gVCFs by default. To skip reblocking, add '\"WholeGenomeReprocessing.WholeGenomeGermlineSingleSample.BamToGvcf.skip_reblocking\": true' to the inputs
* Added WGS plumbing tests for dragen_maximum_quality_mode and dragen_functional_equivalence_mode
* Moved Dragmap docker to WARP and updated to follow repo's best practices
* Added Xmx flag (maximum heap size) to all tasks with java commands
* Added option to allow empty ref_alt file for running BWA mem with masked reference
* Added plumbing input JSON for masked reference
* Updated the SumFloats task used in WholeGenomeGermlineSingleSample.wdl to use python3 instead of python2

# 2.5.0
2021-10-18

* Updated the Whole Genome Germline Single Sample (WGS) workflow to include new DRAGEN-GATK functionality; the default pipeline remains unchanged. Read more about the DRAGEN-GATK mode in the [WGS Overview](https://broadinstitute.github.io/warp/docs/Pipelines/Whole_Genome_Germline_Single_Sample_Pipeline/README)
* Added optional BQSR outputs
* Added a new Docker image for GATK v4.2.2.0 for variant calling in DRAGEN mode

# 2.4.0
2021-10-06

* Updated VerifyBamID to use AppSec base image
* Change GoTC image to Samtools specific image in CramToUnmappedBams and Utilities
* Changed GoTC image to GATK specific image in GermlineVariantDiscovery
* Changed GoTC image to SAMTOOLS/PICARD/BWA specific image in Alignment

# 2.3.9
2021-09-22

* Updated Utilities.wdl task definitions to include a new ErrorWithMessage task that is NOT used in the WholeGenomeReprocessing pipeline.

# 2.3.8
2021-08-02

* Increased the version number to make new release tag for Dockstore

# 2.3.7
2021-06-22

* Removed duplicate MarkDuplicatesSpark task from BamProcessing
* Removed duplicate Docker image from CheckPreValidation task in QC

# 2.3.6
2021-06-01 

* Removed deprecated parameter PAIRED_RUN from MergeBamAlignment

# 2.3.5
2021-03-17

* Promoted VariantCalling to be a top-level workflow

# 2.3.4
2021-02-22

* Added SORTING_COLLECTION_SIZE_RATIO as an optional task input to MarkDuplicates

# 2.3.3
2021-02-08

* Calculate java memory value from the optional memory input value for CramToUnmappedBams java tasks

# 2.3.2
2021-02-02

* Minor changes to support CramToUnmappedBams as an independent versioned pipeline
    * Changed path of the relative import
    * Added 'base_file_name' as an input to CramToUnmappedBams

# 2.3.1
2020-12-21

* Passed an input bam index to several subworkflows, so the pipeline passes on singularity for sharded BQSR

# 2.3.0
2020-12-16

* Fixed error in relative import statement in Alignment subworkflow.
* Fixed syntax bug in Alignment task SamToFastqAndBwaMemAndMba

# 2.2.0
2020-10-20

* Updated GATK docker image for all tasks to [GATK 4.1.8.0](https://github.com/broadinstitute/gatk/releases/tag/4.1.8.0).
    * Numerous bug fixes and improvements
* Updated Picard docker image for all tasks to [2.23.8](https://github.com/broadinstitute/picard/releases/tag/2.23.8).
* Updated samtools to version [1.11](https://github.com/samtools/samtools/releases/tag/1.11).  Primarily for improved compression of cram files.

# 2.1.0
2020-08-18

* Added a meta section with 'allowNestedInputs' set to 'true' to allow the workflow to use with nested inputs with Cromwell 52

# 2.0.2
2020-07-31

* Update various tasks to stop using phusion/baseimage:latest docker image (it has been removed).  Start using a Google-hosted base image in it's stead.

# 2.0.1
2020-07-15

* Remove GetBWAVersion as a task and moved it to SamToFastqAndBwaMemAndMba

# 2.0 
2020-05-13

### Breaking changes to the structure of pipeline inputs. 
* Changes to the inputs included with the dna seq single sample references struct:
    * Removed 'fingerprint_genotypes_file' and 'fingerprint_genotypes_index' from bundle and made these optional pipeline inputs
    * Removed 'haplotype_scatter_count' and 'break_bands_at_multiples_of' from bundle and added these to a separate 'VariantCallingScatterSettings' struct
    * Added 'haplotype_database_file' to the references bundle as a non-optional file
#### Additional changes
* Added ability to convert input_cram to bam using an (optional) alternate reference than the standard.
* Fixed bug in CramToUnmappedBams.RevertSam where it was not reverting the OA tag
* Updated CramToUnmappedBams to properly use the output_map file to support testing.
* Renamed GermlineSingleSampleReferences to DNASeqSingleSampleReferences
* Updated shared tasks to support the new TargetedSomaticSingleSample pipeline

# 1.4
2020-03-05

* Added 'additional_disk' parameter to accommodate larger samples that have steps that run out of disk.

# 1.3
2019-12-03

* Modified the underlying WholeGenomeGermlineSingleSample pipeline to use an up-to-date set of contamination resource files for VerifyBamId.

# 1.2
Adjusted memory parameters to avoid Google's new e2 instances because there are not enough machines to satisfy our production use case.

# 1.1
This update is the result of a a major update to the WholeGenomeGermlineSingleSample pipeline.
We are jumping forward several versions of Picard, from version [2.16.0](https://github.com/broadinstitute/picard/releases/tag/2.16.0) to [2.20.4](https://github.com/broadinstitute/picard/releases/tag/2.20.4)

# 1.0
Initial release of the WholeGenomeReprocessing pipeline
