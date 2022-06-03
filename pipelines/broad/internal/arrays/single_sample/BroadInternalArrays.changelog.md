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