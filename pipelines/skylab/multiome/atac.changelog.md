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

