# 1.0.22
2024-10-28 (Date of Last Commit)

* Updated the docker in the ValidateVCF task; this does not affect this pipeline

# 1.0.21
2024-09-06 (Date of Last Commit)

* Updated to GATK version 4.6.0.0

# 1.0.20
2024-08-02 (Date of Last Commit)

* The ubuntu_16_0_4 docker image version was pinned instead of using the latest tag; this does not affect the outputs of the pipeline

# 1.0.19
2024-06-12 (Date of Last Commit)

* ValidateVcfs is more robust to larger inputs; this does not affect this pipeline

# 1.0.18
2024-07-00 (Date of Last Commit)

* Updated tasks GermlineVariantDiscovery.wdl and QC.wdl to allow multi-cloud dockers; this does not affect this pipeline

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
2023-12-08 (Date of Last Commit)

* ValidateVcfs now has optional memory parameter; this does not affect this pipeline

# 1.0.13
2023-11-29 (Date of Last Commit)

* ValidateVcfs can now take in VCF as calling_interval_list that is in a separate location from its index; this does not affect this pipeline

# 1.0.12
2023-03-30 (Date of Last Commit)
* CheckFingerprint can allow LOD 0

# 1.0.11
2022-11-09 (Date of Last Commit)

* Updated to GATK version 4.3.0.0

# 1.0.10
2022-11-08 (Date of Last Commit)

* Removed task MakeOptionalOutputBam in Utilities.wdl, this update has no effect on this pipeline
* Updated task IngestOutputsToTDR in InternalTasks.wdl with new docker tag to accommodate changes for BroadInternalArrays pipeline. Change has no effect on this pipeline.
* Updated task IngestOutputsToTDR in InternalTasks.wdl with new optional input variable. This update has no effect on this pipeline.
* Updated task FormatArraysOutputs in InternalArrraysTasks.wdl with new docker tag to accommodate changes for BroadInternalArrays pipeline. Change has no effect on this pipeline.
* Removed workspace_bucket workflow parameter from BroadInternalArrays and BroadInternalImputation.

# 1.0.9
2022-09-30 (Date of Last Commit)

* Updated Picard-Python Docker image in Utilities.wdl to fix vulnerabilities.
* Updated task IngestOutputsToTDR with GCR images instead of Dockerhub.

# 1.0.8
2022-09-07 (Date of Last Commit)

* Updated task MakeOptionalOutputBam in Utilities.wdl, this update has no effect on this pipeline
* Updated task IngestOutputsToTDR in InternalTasks.wdl with new docker tag to accommodate changes for BroadInternalArrays pipeline. Change has no effect on this pipeline.

# 1.0.7
2022-06-16 (Date of Last Commit)

* Changed the name of the task QC.CheckFingerprint to QC.CheckFingerprintTask. This prevents a naming conflict in the updated scala tests.

# 1.0.6
2022-06-01 (Date of Last Commit)

* Renamed the CompareVCFs task in VerifyCheckFingerprint.wdl to CompareVcfsAllowingQualityDifferences, this update has no effect on this pipeline

# 1.0.5
2022-05-19 (Date of Last Commit)

* Patch security vulnerability in arrays-picard-private docker image
* Update arrays internal tasks, this update has no effect on this pipeline

# 1.0.4
2022-04-19 (Date of Last Commit)

* Updated to Picard version 2.26.10 and GATK version 4.2.6.1 to address log4j vulnerabilities

# 1.0.3
2022-04-14 (Date of Last Commit)

* Update task CopyFilesFromCloudToCloud in Utilities.wdl, this update has no effect on this pipeline

# 1.0.2
2022-03-24 (Date of Last Commit)

* Update base image for picard-private docker image
* Add gsutil to PATH in picard-private docker image

# 1.0.1
2022-02-01 (Date of Last Commit)

* Fixed bug where fingerprint VCF index file was not passed to the CheckFingerprint task itself
* Addressed memory usage in CheckFingerprint task to allow sufficient headroom for the VM

# 1.0.0
2022-01-14 (Date of Last Commit)

* Initial Release of CheckFingerprint pipeline.
