# 3.0.3
2022-01-20 (Date of Last Commit)

* Increased the disk space in CalibrateDragstrModel task
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
