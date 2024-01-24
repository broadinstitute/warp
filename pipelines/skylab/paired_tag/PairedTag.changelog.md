# 0.0.5
2024-01-24 (Date of Last Commit)

* Added task GetNumSplits before FastqProcess task to determine the number of splits based on the bwa-mem2 machine specs
* Added error message the BWAPairedEndAlignment task to ensure that the number of splits equals to number of ranks

# 0.0.4
2024-01-18 (Date of Last Commit)

* Increased memory for MergeStarOutputs in StarAlign.wdl, RunEmptyDrops in RunEmptyDrops.wdl, OptimusH5ad in H5adUtils.wdl and GeneMetrics in Metrics.wdl
* Added the --soloMultiMappers flag as an optional input to the StarSoloFastq task in the StarAlign.wdl
* Added a check of read2 length and barcode orientation to the demultiplexing step of the pipeline; this task now checks read2 length, performs demultiplexing or trimming if necessary, and checks barcode orientation

# 0.0.3
2024-01-05 (Date of Last Commit)

* Added a new option for the preindex boolean to add cell barcodes and preindex sample barcode to the BB tag of the BAM
* Added new functionality for the ATAC workflow to use BB tag of BAM for SnapATAC2
* Modified the STARsoloFastq task in the StarAlign.wdl so STARsolo can run different types of alignments in a single STARsolo command depending on the counting_mode

# 0.0.2
2023-12-20 (Date of Last Commit)

* Updated the ATAC WDL for the Multiome BWAPairedEndAlignment and MergedBAM tasks

# 0.0.1
2023-12-18 (Date of Last Commit)

* Initial release of the PairedTag pipeline

