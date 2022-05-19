# 1.0

2022-05-19

## Detailed Release Notes

Initial internal release of the single-sample Imputation pipeline's wrapper WDL that executes Imputation.wdl and pushes workflow outputs back to two tables in a TDR dataset - ImputationOutputsTable and ImputationOutputsWideTable.

ImputationOutputsWideTable is an additional table that captures imputed single sample vcfs - one sample per row with chip_well_barcode - unfurled from the `Array[File]` format that is output from the workflow by default. 

* Format outputs of Imputation.wdl and write to a .tsv file.
* Use .tsv file as input to following task that ingests the data into a desginated tables within a TDR dataset (via configured workflow inputs).