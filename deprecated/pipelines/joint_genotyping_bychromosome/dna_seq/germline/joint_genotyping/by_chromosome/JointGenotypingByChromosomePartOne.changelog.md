# 1.5.3
2025-03-06 (Date of Last Commit)

* The JointGenotypingByChromosomePartOne workflow is deprecated and is no longer supported.

# 1.5.2
2024-11-04 (Date of Last Commit)

* Updated to GATK version 4.6.1.0

# 1.5.1
2024-09-10 (Date of Last Commit)

* Update to BGE filtering options in JointCalling; this has no effect on this pipeline

# 1.5.0
2024-09-06 (Date of Last Commit)

* Updated to GATK version 4.6.0.0 
* Updated how no-calls are represented in the output VCFs (0/0 -> ./.) which also changes some annotations in the VCF

# 1.4.12
2023-12-18 (Date of Last Commit)

* Updated to GATK version 4.5.0.0

# 1.4.11
2023-09-08 (Date of Last Commit)

* Added option to hard filter sites outside of provided interval list to HardFilterAndMakeSitesOnlyVcf task

# 1.4.10
2023-06-29 (Date of Last Commit)

* Added extra_args input to the SplitIntervalList task to support the JointGenotypingTasks.wdl

# 1.4.9
2023-05-23 (Date of Last Commit)

* Made disk and memory available as inputs to the JointGenotypingTasks.wdl.

# 1.4.8
2023-04-05 (Date of Last Commit)

* Updated a check for an empty file produced in the GetFingerprintingIntervalIndices task

# 1.4.7
2022-11-04 (Date of Last Commit)

* Updated GATK verison to 4.3.0.0

# 1.4.6
2022-08-26 (Date of Last Commit)

* Added task to JointGenotypingTasks.wdl for UltimaGenomicsJointGenotyping pipeline. This has no effect on this pipeline.

# 1.4.5
2022-06-01 (Date of Last Commit)

* Added inputs to the GenotypeGVCFs task to support the UltimaGenomicsJointGenotyping.wdl

# 1.4.4
2022-04-19 (Date of Last Commit)

* Updated to Picard version 2.26.10 and GATK version 4.2.6.1 to address log4j vulnerabilities

# 1.4.3
2022-04-14 (Date of Last Commit)

* Remove annotationDB files from per chromosome in JointGenotyping

# 1.4.2
2021-11-10

* Task wdls used by JointGenotypingByChromosomePartOne were updated with changes that don't affect JointGenotypingByChromosomePartOne wdl
* Added Xmx flag (maximum heap size) to all tasks with java commands

# 1.4.1
2020-12-16

Modified JointGenotyping 'GatherTranches' task to fix validation error reported by miniwdl (duplicately named input and output)

# 1.4.0
2020-09-25

* Updated GATK docker image for all tasks from [GATK 4.1.1.0](https://github.com/broadinstitute/gatk/releases/tag/4.1.1.0) to [GATK 4.1.8.0](https://github.com/broadinstitute/gatk/releases/tag/4.1.8.0).
    * Numerous bug fixes and improvements, including better VQSR reliability for allele-specific filtering.
    * Addition of new Scattered CrosscheckFingerprints to pipeline

# 1.3.0
2020-09-16

Fix SplitIntervals task so that it works for abutting WGS intervals. Increase by chromosome scatter count appropriately for WGS as fallback.

# 1.2.1
2020-08-28

* Added 'allowNestedInputs: true' metadata parameter to wdl to support Cromwell version 48 and above

# 1.2
Joint Genotyping now accepts a flag to tell VariantRecalibrator to use non-allele-specific annotations. VCFs without allele-specific annotations can now be processed by providing `"JointGenotyping.use_allele_specific_annotations"` in the inputs JSON.

The `"JointGenotypingByChromosomePartOne.GenotypeGVCFs.allow_old_rms_mapping_quality_annotation_data"` input has been added to allow for GVCFs called by the GATK3 HaplotypeCaller to be joint called. GVCFs called with the GATK4 HaplotypeCaller do not need this option set to `true`.

# 1.1
Converging WGS and Exome Joint Genotyping into one all-powerful workflow!

# 1.0
Initial release of the ExomeJointGenotypingByChromosomePartOne pipeline
