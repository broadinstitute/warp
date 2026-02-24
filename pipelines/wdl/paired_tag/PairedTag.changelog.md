# 2.1.10
2026-01-22 (Date of Last Commit)

* Added a new, defaulted input cellbender_memory_GB to Optimus; this does not affect the outputs of this pipeline
* Added a task level input, mem_size, to StarSoloFastq to expose memory settings; this does not affect the outputs of this pipeline

# 2.1.9
2025-09-19 (Date of Last Commit)

* Resolved merge conflicts and reorganize WDL pipelines into unified directory

# 2.1.8
2025-08-15 (Date of Last Commit)

* Added an option to input an aligned ATAC BAM file to the PairedTag pipeline; this allows users to skip the ATAC alignment step if they already have an aligned BAM file

# 2.1.7
2025-08-06 (Date of Last Commit)

* Added MaskPeakCallingMetrics task to remove innappropriate library level peak calling metrics from PairedTag

# 2.1.6
2025-07-31 (Date of Last Commit)

* Added reference genome/GTF headers to fragment file via new string inputs; this change does not affect this pipeline

# 2.1.5
2025-07-16 (Date of Last Commit)

* Added reference genome/GTF headers to fragment file via new string inputs; this change does not affect this pipeline

# 2.1.4
2025-06-18 (Date of Last Commit) 

* Added the exclude_chroms input parameter to the snap.metrics.tsse function in the CreateFragmentFile task of the ATAC pipeline; which is a list of chromosomes to exclude in per cell metric computation. The default value is [chrM, M]

# 2.1.3
2025-06-06 (Date of Last Commit)

* Added mito_list (a list of strings) as an input parameter to CreateFragmentFile task of the ATAC pipeline. This specifies the chromosome names considered mitochondrial DNA. The default value is [chrM, M]
* Removed quotes from bootDiskSizeGb in RunEmptyDrops task to be compatible with Google Batch; this does not affect the outputs of the pipeline
* Increased the ulimit in the following tasks: CalculateCellMetrics, CalculateGeneMetrics CalculateUMIsMetrics; this does not affect the outputs of the pipeline

# 2.1.2
2025-05-27 (Date of Last Commit)

* Increased the ulimit in the STARsoloFastq task in the StarAlign.wdl to 10000; this does not affect the outputs of the pipeline

# 2.1.1

2025-04-15 (Date of Last Commit)

* Refactored peak calling task to be called from PeakCalling.wdl

# 2.1.0

2025-04-02 (Date of Last Commit)
* Refactored the STAR alignment step in Optimus and removed tasks FastqProcessing and MergeSortBamFiles; we are no longer sharding. We are now running one instance of STAR
* Removed MergeStarOutput tasks from Optimus pipeline; added necessary parts of MergeStarOutput task to the STAR alignment step (STARsoloFastq). Additional outputs added to STARsoloFastq task as a result; this includes row_index, col_index, sparse_counts, library_metrics, mtx_files, filtered_mtx_files and cell_reads_out
* Updated the STAR docker image to include Samtools and Python


# 1.10.2 
2025-02-25 (Date of Last Commit)

* Updated the SnapATAC2 docker image to the latest SnapATAC2, allowing for future peak calling implementation
* Updated the warp-tools docker image to include an update to the GroupQCs function in sctools; this does not affect the outputs of the pipeline
* Added reference information to the BAM headers

# 1.10.1
2025-02-03 (Date of Last Commit)

* Added an optional PeakCalling task to the ATAC workflow; this does not affect the outputs of the pipeline
* Added a boolean variable run_peak_calling to the Multiome pipeline; default is false and this does not affect the outputs of the pipeline

# 1.10.0
2025-01-21 (Date of Last Commit)

* Added a boolean variable is_slidetags; default is false, but set to true if Slide-Tags pipeline is calling Optimus
* Added reference_gtf_file to the output h5ad unstructured metadata 
* Added the fragment file CSI index as workflow output
* Updated the default STARsolo multimapping parameter to the EM tehcnique

# 1.9.0
2024-12-05 (Date of Last Commit)

* Added an optional task to the Optimus.wdl that will run CellBender on the Optimus output h5ad file

# 1.8.4
2024-12-3 (Date of Last Commit)

* Fixed a bug in the StarSoloFastq task that caused the pipeline to not output a UniqueAndMult-Uniform.mtx when --soloMultiMappers Uniform was passed to STAR

# 1.8.3
2024-11-22 (Date of Last Commit)

* Added bam validation in the StarSoloFastq task; this does not affect the outputs of the pipeline
* Updated the warp-tools docker; this update changes the way gene_names are identified when creating gene expression h5ad files

# 1.8.2
2024-11-12 (Date of Last Commit)

* Renamed the ATAC workflow library metric percent_target to atac_percent_target for compatibility with downstream tools
* Added more disk and memory to the ParseBarcodes task

# 1.8.1
2024-11-04 (Date of Last Commit)

* Updated the tabix flag in JoinMultiomeBarcodes task in H5adUtils.wdl to use CSI instead of TBI indexing, which supports chromosomes larger than 512 Mbp; this task should not affect the Paired-Tag pipeline

# 1.8.0
2024-10-23 (Date of Last Commit)

* Updated the workflow to include a new expected_cells input parameter describing the number of cells used as input to the library preparation; this is passed to both the ATAC workflows and Optimus workflows and the default is set to 3000 cells
* Updated the ATAC library CSV and the Gene Expression library CSV to be consistent in file naming convention and to have similar case for metric names
* Added a new metric to the ATAC library CSV to calculate percent_target, which is the number of estimated cells by SnapATAC2 divided by expected_cells input
* Updated the ATAC fragment file output so that it is bgzipped
* Updated memory settings for PairedTag Utils

# 1.7.1
2024-10-18 (Date of Last Commit)

* Removed the underscore of the NHashID in the ATAC library metrics CSV

# 1.7.0
2024-09-24 (Date of Last Commit)

* Added a python implementation of DoubletFinder to calculate doublet scores in gene expression data; percent doublets are now available as a library-level metric and individual doublet scores for cell barcodes are in the h5ad
* Updated gene_names in the final h5ad to be unique

# 1.6.1
2024-09-11 (Date of Last Commit)

* Updated warp-tools docker which added create_h5ad_snss2.py to the docker image. This change does not affect the PairedTag pipeline

# 1.6.0
2024-08-02 (Date of Last Commit)

* Updated the SnapATAC2 docker to include v2.7.0; the pipeline will now produce a library-level summary metric CSV for the BAM.

# 1.5.0
2024-08-06 (Date of Last Commit)

* Updated the warp-tools docker to calculate mitochondrial reads from unique reads in cell and gene metrics; these metrics are in the cell and gene metrics CSV as well as h5ad

# 1.4.1
2024-08-02 (Date of Last Commit)

* The ubuntu_16_0_4 docker image version was pinned instead of using the latest tag; this does not affect the outputs of the pipeline

# 1.4.0
2024-07-25 (Date of Last Commit)

* Updated the warp-tools docker image to add TSO metrics to the output h5ad and metric CSV files
* Update the library-level metrics to include new TSO metrics and NHashID descriptor

# 1.3.1
2024-07-18 (Date of Last Commit)

* The atac.wdl was refactored into its own directory under the pipelines/wdl directory; this change does not impact the PairedTag outputs

# 1.3.0
2024-07-11 (Date of Last Commit)

* Updated the PairedTag.wdl to run on Azure. cloud_provider is a new, required input.

# 1.2.0
2024-07-09 (Date of Last Commit)

* Added new optional input parameter of nhash_id, an optional identifier for a library aliquot that is echoed in the  workflow fragment h5ad, the Optimus workflow gene expression h5ad (in the data.uns), and the Optimus gene expression library metrics CSV output; default is set to null
* Added test statements again for GH action (to release from develop). Will probably revert

# 1.1.0
2024-06-28 (Date of Last Commit)

* Updated the STARsolo parameters for estimating cells to Emptydrops_CR
* Added an optional input for expected cells which is used for metric calculation

# 1.0.0
2024-06-26

* Official release of the PairedTag pipeline

# 0.7.0
2024-05-20

* Updated SnapATAC2 docker and tasks to run SnapATAC v2.6.3 
* Added testing infrastructure for paired-tag plumbing data and example data sets


# 0.6.1
2024-05-14 (Date of Last Commit)

* Updated the demultiplex task so that some intermediate input names have been renamed. There is no change to the outputs.

# 0.6.0
2024-05-10 (Date)

* Updated the UPStools docker to version 2.0.0 which reinstates a barcode orientation check script
* Updated the demultiplex task so that output FASTQ files are renamed according to the original input file name, avoiding naming collisions in downstream tasks

# 0.5.1
2024-04-12 (Date of Last Commit)

* Updated the input parameters for STARsolo in STARsoloFastq task. These include the parameters: soloCBmatchWLtype, soloUMIdedup and soloUMIfiltering
* Added "Uniform" as the default string for STARsolo multimapping parameters

# 0.5.0
2024-04-03 (Date of Last Commit)

* Modified adaptor trimming to trim last 3 bp (instead of first) in read2 length is 27 bp and preindex is false

# 0.4.1
2024-03-26 (Date of Last Commit)

* Updated the median umi per cell metric for STARsolo library-level metrics

# 0.4.0
2024-03-15 (Date of Last Commit)

* Added cell metrics to the library-level metrics

* Updated the docker for the MergeStarOutput task to include STARsolo v2.7.11a and custom scripts to create a uniform matrix file and scripts to collect library-level metrics from STARsolo output

* Modified the MergeStarOutput to call a custom script for creating a uniform matrix file (mtx) from individual shard mtx files and to create a filtered matrix from the uniform matrix with STARsolo
# 0.3.0

2024-03-01 (Date of Last Commit)

* Added the gene expression library-level metrics CSV as output of the Paired-tag pipeline; this is produced by the Optimus subworkflow

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

