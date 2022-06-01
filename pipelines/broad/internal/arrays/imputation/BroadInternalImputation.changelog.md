# 1.0.2
2022-06-01 (Date of Last Commit)

* Renamed the CompareVCFs task in VerifyImputation.wdl to CompareVcfsAllowingQualityDifferences, this update has no effect on this pipeline

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