# 1.4.0
2020-09-25

* Updated GATK docker image for all tasks from [GATK 4.1.1.0](https://github.com/broadinstitute/gatk/releases/tag/4.1.1.0) to [GATK 4.1.8.0](https://github.com/broadinstitute/gatk/releases/tag/4.1.8.0).
    * Numerous bug fixes and improvements, including better VQSR reliability for allele-specific filtering.
* Updated picard-cloud docker image from a custom branch to [Picard 2.23.0](https://github.com/broadinstitute/picard/releases/tag/2.23.0)

# 1.3.0
2020-09-16

Fix SplitIntervals task so that it works for abutting WGS intervals. Increase by chromosome scatter count appropriately for WGS as fallback.

# 1.2.1
Added 'allowNestedInputs: true' metadata parameter to wdl to support Cromwell version 48 and above

# 1.2
Joint Genotyping now accepts a flag to tell VariantRecalibrator to use non-allele-specific annotations. VCFs without allele-specific annotations can now be processed by providing `"JointGenotyping.use_allele_specific_annotations"` in the inputs JSON.

# 1.1
Converging WGS and Exome Joint Genotyping into one all-powerful workflow!

# 1.0
Initial release of the ExomeJointGenotypingByChromosomePartTwo pipeline
