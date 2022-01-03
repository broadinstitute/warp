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
