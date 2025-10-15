# 1.2.1
2025-10-09 (Date of Last Commit)

* Modified the ReblockGVCF.wdl to use bash's basename instead of WDL's basename; this does not affect the outputs of this pipeline


# 1.2.0
2025-03-17 (Date of Last Commit)

* Updated Picard docker for CompareMetrics from 2.20.4-SNAPSHOT to 2.26.10-SNAPSHOT to fix security vulnerability. The outputs of the CollectAggregationMetrics task have changed slightly: .alignment_summary_metrics now have additional columns SD_READ_LENGTH, MEDIAN_READ_LENGTH, MAD_READ_LENGTH, MIN_READ_LENGTH, MAX_READ_LENGTREADS_ALIGNED_IN_PAIRS, PCT_SOFTCLIP, PCT_HARDCLIP, AVG_POS_3PRIME_SOFTCLIP_LENGTH. Additionally, a histogram is now included in the .alignment_summary_metrics file

# 1.1.4
2025-02-21 (Date of Last Commit)

* Updated HaplotypeCaller_GATK4_VCF to use MEM_SIZE and MEM_UNIT; this does not affect the outputs of this pipeline

# 1.1.3
2024-12-05 (Date of Last Commit)

* Updated the name of the output for ReblockGVCFs; this does not affect this pipeline

# 1.1.2
2024-11-04 (Date of Last Commit)

* Updated to GATK version 4.6.1.0

# 1.1.1
2024-10-28 (Date of Last Commit)

* Updated the docker in the ValidateVCF task; this does not affect this pipeline

# 1.1.0
2024-09-06 (Date of Last Commit)

* Updated to GATK version 4.6.0.0
* Minor changes to the output VCF at some sites due to the changes in HaplotyperCaller

# 1.0.20
2024-08-02 (Date of Last Commit)

* The ubuntu_16_0_4 docker image version was pinned instead of using the latest tag; this does not affect the outputs of the pipeline

# 1.0.19
2024-06-12 (Date of Last Commit)

* ValidateVcfs is more robust to larger inputs; this does not affect this pipeline

# 1.0.18
2024-07-09 (Date of Last Commit)

* Updated GermlineVariantDiscovery, BamProcessing, DragenTasks, Qc, and Utilities tasks to allow multi-cloud dockers. This change does not affect this pipeline

# 1.0.17
2024-07-01 (Date of Last Commit)

* CalculateReadGroupChecksum requires more memory and disk; this does not affect this pipeline

# 1.0.16
2024-03-26 (Date of Last Commit)

* ValidateVcfs requires less memory when run without interval list; this does not affect this pipeline

# 1.0.15
2023-12-18 (Date of Last Commit)

* Updated to GATK version 4.5.0.0.

# 1.0.14
2023-12-14 (Date of Last Commit)

* Updated GATK for Reblock task to version 4.5.0.0
* Added options to Reblock task to remove annotations and move filters to genotype level

# 1.0.13
2023-12-08 (Date of Last Commit)

* ValidateVcfs now has optional memory parameter; this does not affect this pipeline

# 1.0.12
2023-11-29 (Date of Last Commit)

* ValidateVcfs can now take in VCF as calling_interval_list that is in a separate location from its index; this does not affect this pipeline

# 1.0.11
2023-11-21 (Date of Last Commit)

* Fixed bug in ReblockGVCFs; this does not affect this pipeline.
* Reverted the VerifyBamID docker image back to the 1.0.9 UltimaGenomicsWholeGenomeGermline pipeline version

# 1.0.10
2023-09-18 (Date of Last Commit)

* ReblockGVCFs can now take in GVCFs that are not in the same location as their index file, this update has no effect on this pipeline.

# 1.0.9
2023-08-16 (Date of Last Commit)

* Updated VerifyBamID docker image in UltimaGenomicsWholeGenomeGermlineTasks.wdl to fix security vulnerabilities, this update has no effect on this pipeline.

# 1.0.8
* CheckFingerprint can allow LOD 0

# 1.0.7
2023-01-13 (Date of Last Commit)

* Updated remaining usese of GATK to verison 4.3.0.0.


# 1.0.6
2022-11-04 (Date of Last Commit)

* Updated GATK verison to 4.3.0.0. This includes several bug fixes in HaplotypeCaller.
* Added UG High Quality interval list to single sample random forest filtering
* Changed cutoff QUAL value in GVCF->VCF conversion from 0 to 30. This removes lower quality sites from single sample random forest filtering and from the final VCF since sites that aren't scored by the random forest are removed in ReblockGVCFs.

# 1.0.5
2022-11-08 (Date of Last Commit)

* Removed task MakeOptionalOutputBam in Utilities.wdl, this update makes this pipeline more robust for large samples
* Updated task IngestOutputsToTDR in InternalTasks.wdl with new docker tag to accommodate changes for BroadInternalArrays pipeline. Change has no effect on this pipeline.
* Updated task IngestOutputsToTDR in InternalTasks.wdl with new optional input variable. This update has no effect on this pipeline.
* Updated task FormatArraysOutputs in InternalArrraysTasks.wdl with new docker tag to accommodate changes for BroadInternalArrays pipeline. Change has no effect on this pipeline.
* Removed workspace_bucket workflow parameter from BroadInternalArrays and BroadInternalImputation.
 
# 1.0.4
2022-09-30 (Date of Last Commit)

* Updated Picard-Python Docker image in Utilities.wdl to fix vulnerabilities.
* Updated task IngestOutputsToTDR with GCR images instead of Dockerhub.

# 1.0.3
2022-09-20 (Date of Last Commit)

* Removed /cromwell_root/ prefix for output file paths in FilterVCF and TrainModel tasks.

# 1.0.2
2022-09-07 (Date of Last Commit)

* Increased disk space in the MakeOptionalOutputBam task in Utilities.wdl
* Updated task IngestOutputsToTDR in InternalTasks.wdl with new docker tag to accommodate changes for BroadInternalArrays pipeline. Change has no effect on this pipeline.

# 1.0.1
2022-06-21 (Date of Last Commit)

* Changed QC.CheckFingerprint to QC.CheckFingerprintTask to avoid a naming conflict in the update scala tests, no effect on this pipeline

# 1.0.0
2022-05-05 (Date of Last Commit)

* Initial Release of UltimaGenomicsWholeGenomeGermline pipeline.
* The UltimaGenomicsWholeGenomeGermline pipeline is an open-source, cloud-optimized workflow created for processing Ultima Genomics Whole Genome Sequenced Germline samples. Overall, the workflow aligns reads to the genome, marks duplicates, calls variants, and calculates quality metrics to produce a CRAM, CRAI, GVCF, filtered VCF, and quality metrics.