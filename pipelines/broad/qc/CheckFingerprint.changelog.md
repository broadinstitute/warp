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
