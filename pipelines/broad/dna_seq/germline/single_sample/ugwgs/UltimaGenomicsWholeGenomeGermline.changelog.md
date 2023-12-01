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