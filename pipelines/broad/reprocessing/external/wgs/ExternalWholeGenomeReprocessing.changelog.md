# 2.0.4
2022-02-02 (Date of Last Commit)

* Changed dragmap base image from Centos to RockyLinux to comply with trivy scans

# 2.0.3
2022-02-01 (Date of Last Commit)

* Increased the disk space in Reblock task
* Increased the disk space in CalibrateDragstrModel task
* Addressed memory usage in CheckFingerprint task to allow sufficient headroom for the VM

# 2.0.2
2022-01-14 (Date of Last Commit)

* Increased the disk space in CalibrateDragstrModel task

# 2.0.1
2021-12-09
* Updated the base image for the Dragmap docker image
* Updated broken dependency in VerifyBamID docker image

# 2.0.0
2021-11-15

* Added an optional step to reblock gVCFs, this step is included by default
    * The ExternalWholeGenomeReprocessing pipeline now outputs reblocked gVCFs by default. To skip reblocking, add '\"ExternalWholeGenomeReprocessing.WholeGenomeReprocessing.WholeGenomeGermlineSingleSample.BamToGvcf.skip_reblocking\": true' to the inputs
* Added WGS plumbing tests for dragen_maximum_quality_mode and dragen_functional_equivalence_mode
* Moved Dragmap docker to WARP and updated to follow repo's best practices
* Added Xmx flag (maximum heap size) to all tasks with java commands
* Added option to allow empty ref_alt file for running BWA mem with masked reference
* Added plumbing input JSON for masked reference
* Updated the SumFloats task used in WholeGenomeGermlineSingleSample.wdl to use python3 instead of python2

# 1.5.0
2021-10-18

* Updated the Whole Genome Germline Single Sample (WGS) workflow to include new DRAGEN-GATK functionality; the default pipeline remains unchanged. Read more about the DRAGEN-GATK mode in the [WGS Overview](https://broadinstitute.github.io/warp/docs/Pipelines/Whole_Genome_Germline_Single_Sample_Pipeline/README)
* Added optional BQSR outputs
* Added a new Docker image for GATK v4.2.2.0 for variant calling in DRAGEN mode

# 1.4.0
2021-10-06

* Updated VerifyBamID to use AppSec base image
* Change GoTC image to Samtools specific image in CramToUnmappedBams and Utilities
* Changed GoTC image to GATK specific image in GermlineVariantDiscovery

# 1.3.9
2021-09-22

* Updated Utilities.wdl task definitions to include a new ErrorWithMessage task that is NOT used in the ExternalWholeGenomeReprocessing pipeline.

# 1.3.8
2021-08-02

* Increased the version number to make new release tag for Dockstore 

# 1.3.7
2021-06-22

* Removed duplicate MarkDuplicatesSpark task from BamProcessing
* Removed duplicate Docker image from CheckPreValidation task in QC

# 1.3.6
2021-06-01 

* Removed deprecated parameter PAIRED_RUN from MergeBamAlignment

# 1.3.5
2021-03-17

* Promoted VariantCalling to be a top-level workflow

# 1.3.4
2021-02-22

* Added SORTING_COLLECTION_SIZE_RATIO as an optional task input to MarkDuplicates

# 1.3.3
2021-02-08

* Calculate java memory value from the optional memory input value for CramToUnmappedBams java tasks

# 1.3.2
2021-02-02

* Minor changes to support CramToUnmappedBams as an independent versioned pipeline
    * Changed path of the relative import
    * Added 'base_file_name' as an input to CramToUnmappedBams

# 1.3.1
2020-12-21

* Passed an input bam index to several subworkflows, so the pipeline passes on singularity for sharded BQSR

# 1.3.0
2020-12-16

* Fixed error in relative import statement in Alignment subworkflow.
* Fixed syntax bug in Alignment task SamToFastqAndBwaMemAndMba

# 1.2.0
2020-10-20

* Updated GATK docker image for all tasks to [GATK 4.1.8.0](https://github.com/broadinstitute/gatk/releases/tag/4.1.8.0).
    * Numerous bug fixes and improvements
* Updated Picard docker image for all tasks to [2.23.8](https://github.com/broadinstitute/picard/releases/tag/2.23.8).
* Updated samtools to version [1.11](https://github.com/samtools/samtools/releases/tag/1.11).  Primarily for improved compression of cram files.

# 1.1.1
2020-10-01

* Removed extra trailing slash in ouput directory from cloud to cloud copy job

# 1.1.0
2020-08-18

* Added a meta section with 'allowNestedInputs' set to 'true' to allow the workflow to use with nested inputs with Cromwell 52

# 1.0.2
2020-07-31

* Update various tasks to stop using phusion/baseimage:latest docker image (it has been removed).  Start using a Google-hosted base image in it's stead.

# 1.0.1
2020-07-15

* Remove GetBWAVersion as a task and moved it to SamToFastqAndBwaMemAndMba

# 1.0
2020-04-30

* Initial release of the ExternalWholeGenomeReprocessing pipeline
