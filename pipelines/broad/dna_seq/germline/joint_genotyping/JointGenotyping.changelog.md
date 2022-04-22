# 1.6.0
2022-04-14 (Date of Last Commit)

* Updated to Picard version 2.26.10 and GATK version 4.2.6.1 to address log4j vulnerabilities
    * ExcessHet values change from removing mid-p correction; now PASS is zero 
    * AS_\*RankSum values change after a fix for Histogram::median() 
    * AS_MQ, AS_FS, AS_SOR values change in some cases from GenomicsDB non-ref assignment fix
    * \* allele no longer gets annotation values
    * Gnarly AN differences: Some ./. go to \*/\*, but those agree with upstream deletions

# 1.5.3
2022-04-12

* Remove annotationDB files from per chromosome in JointGenotyping

# 1.5.2
2021-11-10

* Updated GenotypeGVCFs to support reblocked GVCFs as inputs
* Added Xmx flag (maximum heap size) to all tasks with java commands

# 1.5.1
2020-12-16

Modified JointGenotyping 'GatherTranches' task to fix validation error reported by miniwdl (duplicately named input and output)

# 1.5.0
2020-10-20

Updated the sample map for WGS JointGenotyping to one with new samples. Turned on validation of NA12878 for both Exome and WGS JointGenotyping.

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

The `"JointGenotyping.GenotypeGVCFs.allow_old_rms_mapping_quality_annotation_data"` input has been added to allow for GVCFs called by the GATK3 HaplotypeCaller to be joint called. GVCFs called with the GATK4 HaplotypeCaller do not need this option set to `true`.

# 1.1
Converging WGS and Exome Joint Genotyping into one all-powerful workflow!

# 1.0
Initial release of the ExomeJointGenotyping pipeline