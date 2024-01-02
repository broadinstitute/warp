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

