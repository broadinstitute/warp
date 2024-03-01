# 0.2.0
2024-02-29 (Date of Last Commit)
* Added mem and disk to inputs of Join Barcodes task of Multiome workflow; does not impact the Paired-tag workflow

# 0.1.0
2024-02-22 (Date of Last Commit)

* Updated StarAlign output metrics to include shard ids, which is called by Optimus
* Remove ref_genome_fasta from Optimus input

# 0.0.7
2024-02-07 (Date of Last Commit)

* Updated the Metrics tasks to exclude mitochondrial genes from reads_mapped_uniquely, reads_mapped_multiple and reads_mapped_exonic, reads_mapped_exonic_as and reads_mapped_intergenic

# 0.0.6
2024-02-01 (Date of Last Commit)

* Add new paired-tag task to parse sample barcodes from cell barcodes when preindexing is set to true

# 0.0.5
2024-01-30 (Date of Last Commit)

* Added task GetNumSplits before FastqProcess ATAC task to determine the number of splits based on the bwa-mem2 machine specs
* Added an error message to the BWAPairedEndAlignment ATAC task to ensure that the number of splits equals the number of ranks
* Added an error message to the BWAPairedEndAlignment ATAC task to ensure that the number of R1s equals to the number of R3s

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

