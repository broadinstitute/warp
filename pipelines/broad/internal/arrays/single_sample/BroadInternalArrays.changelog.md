# 1.1.2
2023-01-13 (Date of Last Commit)

* Updated remaining uses of GATK to version 4.3.0.0

# 1.1.1
2022-12-16 (Date of Last Commit)

* Updated to GATK version 4.3.0.0

# 1.1.0
2022-11-21 (Date of Last Commit)

* add `lab_batch` as required input paramater to workflow BroadInternalArrays
* add `lab_batch` to task FormatArraysOutputs

# 1.0.9 
2022-11-08 (Date of Last Commit)

* remove workspace_bucket parameter from workflow inputs and IngestOutputsToTDR

# 1.0.8
2022-10-12 (Date of Last Commit)

* Removed task MakeOptionalOutputBam in Utilities.wdl, this update has no effect on this pipeline
* Added `in_load_tag` as optional input to task IngestOutputsToTDR. This update has no effect on this pipeline.

* New GCR image tag in task FormatArraysOutputs and task IngestOutputsToTDR.

# 1.0.7
2022-09-30 (Date of Last Commit)

* Updated Picard-Python Docker image in Utilities.wdl to fix vulnerabilities.
* Updated tasks FormatArraysOutputs and IngestOutputsToTDR with GCR images instead of Dockerhub.

# 1.0.6
2022-09-20 (Date of Last Commit)

* Updated call to IngestOutputsToTDR to remove 'prefix_column'. Python script has been updated and not longer requires this input parameter.
* Update task IngestOutputsToTDR to not require 'prefix_column'. Python script has been updated and not longer requires this input parameter.

* Update task FormatArraysOutputs with new docker tag.
* Update task IngestOutputsToTDR with new docker tag.
* Update tasks FormatArraysOutputs and IngestOutputsToTDR with GCR image instead of DockerHub image.

# 1.0.5
2022-08-29 (Date of Last Commit)

* Updated call to IngestOutputsToTDR to pass in column names to be used for user action in command block. Python script in task was updated to a new version containing a new required command line parameter, 'prefix_column'

# 1.0.4
2022-07-15 (Date of Last Commit)

* Updated task MakeOptionalOutputBam in Utilities.wdl, this update has no effect on this pipeline

# 1.0.3
2022-06-21 (Date of Last Commit)

* Changed QC.CheckFingerprint to QC.CheckFingerprintTask to avoid a naming conflict in the update scala tests, no effect on this pipeline

# 1.0.2
2022-06-03 (Date of Last Commit)

* Renamed the CompareVCFs task in VerifyIlluminaGenotypingArray.wdl to CompareVcfsAllowingQualityDifferences, this update has no effect on this pipeline

# 1.0.1
2022-06-02 (Date of Last Commit)

* Fix level of relative import of a pipeline WDL to allow for warp testing framework to operate properly

# 1.0.0
2022-05-19 (Date of Last Commit)

## Initial pipeline release

* Initial internal release of the single-sample Arrays pipeline's wrapper WDL that executes Arrays.wdl and pushes workflow outputs back to a table in a TDR dataset - ArraysOutputsTable.
* Format outputs of Arrays.wdl and write to a .tsv file.
* Use .tsv file as input to following task that ingests the data into a desginated table within a TDR dataset (via configured workflow inputs).