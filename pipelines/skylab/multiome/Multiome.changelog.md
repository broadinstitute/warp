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

