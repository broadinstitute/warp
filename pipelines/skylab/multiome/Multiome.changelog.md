# 3.0.1
2023-12-12 (Date of Last Commit)
* Downgraded Cell Bender from v0.3.1 to v0.3.0
  
# 3.0.0
2023-11-22 (Date of Last Commit)

* Added gene expression barcodes to the Multiome ATAC fragment file
* Updated the JoinBarcodes task to bgzip and tabix the final ATAC fragment file
* Added the tabix index file as an output to Multiome

# 2.3.3
2023-11-21 (Date of Last Commit)

* Added the latest warp-tools docker to tasks in the Metrics, FastqProcessing and H5adUtils wdls; this incorporates new input parameter for number of output fastq files to fastqprocess

# 2.3.2
2023-11-20 (Date of Last Commit)
* Added an optional task to the Multiome.wdl that will run CellBender on the Optimus output h5ad file

# 2.3.1
2023-11-20 (Date of Last Commit)

* Added the latest warp-tools docker to the Metrics task; this allows use of REFSEQ references

# 2.3.0
2023-11-03 (Date of Last Commit)

* Updated the Metrics task so that Cell Metrics and Gene Metrics now calculate intronic, intronic_as, exonic, exonic_as, and intergenic metrics from unique reads only using the NH:i:1 tag in the BAM

# 2.2.2
2023-10-20 (Date of Last Commit)

* Removed Dropna from JoinBarcodes subtask of the H5adUtils task, which was causing the JoinBarcodes to fail for some gene expression matrices

* Updated path to Multiome whitelists to reflect location in public storage.

# 2.2.0
2023-10-05 (Date of Last Commit)
* Added a JoinMultiomeBarcodes task to the H5adUtils that adds a column in the ATAC and Optimus output h5ad linking gene expression and ATAC barcodes

# 2.1.0
2023-09-21 (Date of Last Commit)
* Added dynamic barcode orientation selection to the ATAC workflow FastqProcess task

# 2.0.0
2023-09-05 (Date of Last Commit)

* Updated Optimus pipeline to include STARsolo v2.7.11a
* Added sF tag to STARsolo aligner parameters
* Updated TagSort tool for Optimus Metrics task to calculate metrics based on the sF tag
* Modified H5adUtils task to include new metrics in the final Optimus h5ad
* Removed the Dropseq metrics task
* Updated the ATAC wdl to optimize the MakeFragmentFile task for cpu, memory, disk, and cpu platform
 
# 1.0.1 
2023-07-23 (Date of Last Commit)

* Added STARsolo v2.7.10b metric outputs as an optional pipeline output and an output of the STARalign and MergeSTAR tasks

* Updated the CountAlignments task in the FeatureCounts.wdl to use a new docker image. This change does not affect the Multiome pipeline

# 1.0.0
2023-06-22 (Date of Last Commit)

* Initial release of the multiome pipeline

