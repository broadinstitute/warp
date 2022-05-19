# 1.0

2022-05-19

## Detailed Release Notes

Initial internal release of the single-sample Arrays pipeline's wrapper WDL that executes Arrays.wdl and pushes workflow outputs back to a table in a TDR dataset - ArraysOutputsTable.

* Format outputs of Arrays.wdl and write to a .tsv file.
* Use .tsv file as input to following task that ingests the data into a desginated table within a TDR dataset (via configured workflow inputs).