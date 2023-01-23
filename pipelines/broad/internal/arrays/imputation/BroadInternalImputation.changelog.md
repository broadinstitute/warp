# 1.1.6
2023-01-13 (Date of Last Commit)

* Updated to GATK version 4.3.0.0

# 1.1.5
2022-12-01 (Date of Last Commit)

* Updated BCFTools/VCFTools Docker image

# 1.1.4
2022-11-09 (Date of Last Commit)

* Updated to GATK version 4.3.0.0
* Remove workspace_bucket parameter from workflow inputs and IngestOutputsToTDR


# 1.1.3
2022-11-10 (Date of Last Commit)

* Added meta section to allowNestedInputs, this will allow task level inputs to be set in the inputs json
* Added `in_load_tag` as optional input to task IngestOutputsToTDR.
* Pass `in_load_tag` in the call to IngestToImputationWideOutputsTable (alias for IngestOutputsToTDR).

* New GCR image tag in task FormatImputationOutputs and task FormatImputationWideOutputs.

# 1.1.2
2022-09-27 (Date of Last Commit)

* Removed task MakeOptionalOutputBam in Utilities.wdl, this update has no effect on this pipeline

# 1.1.1
2022-09-30 (Date of Last Commit)

* Updated Picard-Python, BCFTools/VCFTools, and Minimac4 Docker images to fix vulnerabilities.
* Updated tasks FormatImputationOutputs, FormatImputationWideOutputs, and IngestOutputsToTDR with GCR images instead of Dockerhub.

# 1.1.0
2022-09-20 (Date of Last Commit)

* Updated call to IngestOutputsToTDR to remove 'prefix_column'. Python script has been updated and not longer requires this input parameter.
* Update task IngestOutputsToTDR to not require 'prefix_column'. Python script has been updated and not longer requires this input parameter.

* Update task FormatImputationOutputs with new docker tag.
* Update task FormatImputationWideOutputs with new docker tag.
* Update task IngestOutputsToTDR with new docker tag.
* Update tasks FormatImputationOutputs, FormatImputationWideOutputs, and IngestOutputsToTDR with GCR image instead of DockerHub image.

# 1.0.9
2022-08-29 (Date of Last Commit)

* Updated call to IngestOutputsToTDR to pass in column names to be used for user action in command block. Python script in task was updated to a new version containing a new required command line parameter, 'prefix_column'

# 1.0.8
2022-08-23 (Date of Last Commit)

* Updated BCFTools/VCFTools docker image

# 1.0.7
2022-08-03 (Date of Last Commit)

* Updated BCFTools/VCFTools Minimac4 Docker images

# 1.0.6
2022-07-15 (Date of Last Commit)

* Updated task MakeOptionalOutputBam in Utilities.wdl, this update has no effect on this pipeline

# 1.0.5
2022-06-10 (Date of Last Commit)

* Add run_task variable as input to TriggerPrsWithImputationTsv task to ensure that task runs after IngestToImputationWideOutputsTable. Without completion of IngestToImputationWideOutputsTable, the CF that is triggered in TriggerPrsWithImputationTsv will fail - adding the additional variable forces the task to wait for the upstream ingest task to complete and generate its output file, which is required as input to the downstream task.

# 1.0.4
2022-06-06 (Date of Last Commit)

* Add timestamp variable to workflow inputs
* Update TriggerPrsWithImputationTsv task with destination filename prefixed with timestamp from workflow input

 # 1.0.3
2022-06-03 (Date of Last Commit)

* Renamed the CompareVCFs task in VerifyImputation.wdl to CompareVcfsAllowingQualityDifferences, this update has no effect on this pipeline

# 1.0.2
2022-06-02 (Date of Last Commit)

* Fix level of relative import of a pipeline WDL to allow for warp testing framework to operate properly

# 1.0.1
2022-05-23 (Date of Last Commit)

* Task (and call) added to copy the outputs of FormatImputationOutputs to an external GCS location which acts as a trigger for a CF (cloud function) for a downstream eMerge process - set up for PRS scoring WDL.


# 1.0.0
2022-05-19 (Date of Last Commit)

## Initial pipeline release

* Initial internal release of the single-sample Imputation pipeline's wrapper WDL that executes Imputation.wdl and pushes workflow outputs back to two tables in a TDR dataset - ImputationOutputsTable and ImputationOutputsWideTable.

* ImputationOutputsWideTable is an additional table that captures imputed single sample vcfs - one sample per row with chip_well_barcode - unfurled from the `Array[File]` format that is output from the workflow by default. 

* Format outputs of Imputation.wdl and write to a .tsv file.
* Use .tsv file as input to following task that ingests the data into a desginated tables within a TDR dataset (via configured workflow inputs).