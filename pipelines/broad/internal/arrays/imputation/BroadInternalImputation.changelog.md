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