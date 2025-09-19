# 2.9.3
2025-09-19 (Date of Last Commit)

* Resolved merge conflicts and reorganize WDL pipelines into unified directory

# 2.9.2
2025-08-15 (Date of Last Commit)

* Added an option to input an aligned ATAC BAM file; this allows users to skip the ATAC alignment step if they already have an aligned BAM file

# 2.9.1
2025-08-06 (Date of Last Commit) 

* Added MaskPeakCallingMetrics task to remove innappropriate peak calling metrics from PairedTag

# 2.9.0
2025-06-18 (Date of Last Commit) 

* Added the exclude_chroms input parameter to the snap.metrics.tsse function in the CreateFragmentFile task; which is a list of chromosomes to exclude in per cell metric computation. The default value is [chrM, M]

# 2.8.0
2025-06-06 (Date of Last Commit) 

* Added mito_list (a list of strings) as an input parameter to CreateFragmentFile task. This specifies the chromosome names considered mitochondrial DNA. The default value is [chrM, M]


# 2.7.2
2025-04-15 (Date of Last Commit) 

* Refactored peak calling task to be called from PeakCalling.wdl

# 2.7.1
2025-02-25 (Date of Last Commit)

* Added a new warning for peak calling step if the probability_threshold is too low, resutling in a null matrix after doublet filtering
* Updated the probability threshold default to 0.5
* Updated the warp-tools docker image to include an update to the GroupQCs function in sctools; this does not affect the outputs of the pipeline
* Added reference information to the BAM header

# 2.7.0
2025-02-03 (Date of Last Commit)

* Added an optional PeakCalling task 
* Added a boolean variable peak_calling; default is false 

# 2.6.0
2025-01-21 (Date of Last Commit)

* Added reference_gtf_file to the output h5ad unstructured metadata
* Added the fragment file CSI index as a workflow output

# 2.5.3
2024-11-22 (Date of Last Commit)

* Updated the warp-tools docker; this update changes the way gene_names are identified when creating gene expression h5ad files; does not impact ATAC workflow

# 2.5.2
2024-11-12 (Date of Last Commit)

* Added memory and disk updates to Multiome JoinBarcodes; this does not impact the ATAC workflow

# 2.5.1
2024-11-12 (Date of Last Commit)

* Renamed the ATAC workflow library metric percent_target to atac_percent_target for compatibility with downstream tools

# 2.5.0
2024-10-23 (Date of Last Commit)

* Updated the tabix flag in CreateFragmentFile task to use CSI instead of TBI indexing, which supports chromosomes larger than 512 Mbp
* Renamed the ATAC workflow library metric percent_target to atac_percent_target for compatibility with downstream tools

# 2.4.0
2024-10-23 (Date of Last Commit)

* Added a new input parameter for atac_expected_cells, which describes the numnber of cells used for the library preparation
* Updated the ATAC library CSV to be consistent in file naming convention and to have similar case for metric names to the Optimus workflow library CSV
* Added a new metric to the ATAC library CSV to calculate percent_target, which is the number of estimated cells by SnapATAC2 divided by expected_cells input
* Updated the ATAC workflow so that the output fragment file is bgzipped by default
* Updated memory settings for PairedTag; does not impact the ATAC workflow


# 2.3.2
2024-10-18 (Date of Last Commit)

* Removed the underscore of the NHashID in the ATAC library metrics CSV

# 2.3.1
2024-09-11 (Date of Last Commit)

* Updated warp-tools docker which added create_h5ad_snss2.py to the docker image. This change does not affect the atac pipeline

# 2.3.0
2024-08-29 (Date of Last Commit)

* Updated the SnapATAC2 docker to include v2.7.0; the pipeline will now produce a library-level summary metric CSV for the BAM. 

* Updated the memory for the CreateFragmentFile task

# 2.2.3
2024-08-02 (Date of Last Commit)

* Updated the warp-tools docker which now includes new metric calculations for mitochondria reads; this does not impact the ATAC workflow

# 2.2.2
2024-08-02 (Date of Last Commit)

* The ubuntu_16_0_4 docker image version was pinned instead of using the latest tag; this does not affect the outputs of the pipeline

# 2.2.1
2024-07-25 (Dat of Last Commit)

* Updated the warp-tools docker image to add TSO metrics to the output h5ad and metric CSV files; this does not impact the ATAC workflow

# 2.2.0
2024-07-11 (Date of Last Commit)

* Updated the atac.wdl to run on Azure. cloud_provider is a new, required input.

# 2.1.0
2024-07-09 (Date of Last Commit)

* Added new optional input parameter of atac_nhash_id, an identifier for a library aliquot that is echoed in the atac fragment metrics h5ad (in the data.uns); default is set to null 
* Added test statements again for GH action (to release from develop). Will probably revert

# 2.0.0
2024-05-20 (Date of Last Commit)

* Updated SnapATAC2 docker to SnapATAC2 v2.6.3; this impacts the workflow output metrics

# 1.2.3
2024-05-14 (Date of Last Commit)

* Updated the demultiplex task so that some intermediate input names have been renamed, this does not impact the ATAC workflow

# 1.2.2
2024-05-10 (Date of Last Commit)

* Updated the Paired-tag Demultiplex task; this does not impact the ATAC workflow

# 1.2.1
2024-04-03 (Date of Last Commit)
* Modified adaptor trimming in Paired-tag WDL; this does not impact ATAC

# 1.2.0
2024-03-15 (Date of Last Commit)

* Added cell metrics to the library-level metrics

* Updated the docker for the MergeStarOutput task to include STARsolo v2.7.11a and custom scripts to create a uniform matrix file and scripts to collect library-level metrics from STARsolo output

* Modified the MergeStarOutput to call a custom script for creating a uniform matrix file (mtx) from individual shard mtx files and to create a filtered matrix from the uniform matrix with STARsolo

# 1.1.8
2024-02-07 (Date of Last Commit)

* Updated the Metrics tasks to exclude mitochondrial genes from reads_mapped_uniquely, reads_mapped_multiple and reads_mapped_exonic, reads_mapped_exonic_as and reads_mapped_intergenic

# 1.1.7
2024-02-01 (Date of Last Commit)

* Add new paired-tag task to parse sample barcodes from cell barcodes when preindexing is set to true; this does not affect the ATAC workflow

# 1.1.6
2024-01-24 (Date of Last Commit)

* Added task GetNumSplits before FastqProcess task to determine the number of splits based on the bwa-mem2 machine specs
* Added an error message to the BWAPairedEndAlignment task to ensure that the number of splits equals the number of ranks
* Added an error message to the BWAPairedEndAlignment task to ensure that the number of R1s equals to the number of R3s

# 1.1.5 
2024-01-10 (Date of Last Commit)

* Added a check of read2 length to the paired-tag pipeline; this does not affect the ATAC workfow

# 1.1.4
2024-01-02 (Date of Last Commit)

* Added functionality for using the ATAC pipeline with paired-tag data, including the option for SnapATAC task to pull cell barcodes from the BB tag of the BAM

# 1.1.3
2023-12-17 (Date of Last Commit)

* Added updated docker to BWAPairedEndAlignment ATAC task to use updated code for distributed bwa-mem2 from Intel
* Removed MergedBAM ATAC and moved BWAPairedEndAlignment ATAC outside of the for loop
* Changed CPU platform to Ice Lake for BWAPairedEndAlignment ATAC task
* Added input parameter input_output_parameter to the Multiome ATAC wdl
* Increased memory for JoinMultiomeBarcodes in H5adUtils 

# 1.1.2
2023-11-21 (Date of Last Commit)

Added the latest warp-tools docker to tasks in the Metrics, FastqProcessing and H5adUtils wdls; this incorporates new input parameter for number of output fastq files to fastqprocess

# 1.1.1
2023-10-20 (Date of Last Commit)
* Removed the dropna from the JoinBarcodes subtask of the H5adUtils WDL; this change does not impact ATAC outputs
* Added dynamic barcode selection to the ATAC FastqProcess task

# 1.0.1
2023-09-05

* Optimized the MakeFragmentFile task for memory, disk, cpu and cpu platform

# 1.0.0
2023-06-22 (Date of Last Commit)

* Initial release of the ATAC pipeline

