# 3.3.0
2024-02-28 (Date of Last Commit)

* Added the gene expression library-level metrics CSV as output of the Multiome pipeline; this is produced by the Optimus subworkflow

# 3.2.0
2024-02-22 (Date of Last Commit)

* Updated StarAlign.MergeStarOutput to add a shard number to the metrics files
* Removed ref_genome_fasta input from Multiome WDL and JSON

# 3.1.3
2024-02-07 (Date of Last Commit)

* Updated the Metrics tasks to exclude mitochondrial genes from reads_mapped_uniquely, reads_mapped_multiple and reads_mapped_exonic, reads_mapped_exonic_as and reads_mapped_intergenic

# 3.1.2
2024-02-01 (Date of Last Commit)

* Add new paired-tag task to parse sample barcodes from cell barcodes when preindexing is set to true; this does not affect the Multiome pipeline


# 3.1.1 
2024-01-30 (Date of Last Commit)

* Added task GetNumSplits before FastqProcess ATAC task to determine the number of splits based on the bwa-mem2 machine specs
* Added an error message to the BWAPairedEndAlignment ATAC task to ensure that the number of splits equal the number of ranks
* Added an error message to the BWAPairedEndAlignment ATAC task to ensure that the number of R1s equals to the number of R3s

# 3.1.0
2024-01-24 (Date of Last Commit)
* Promote aligner_metrics from Optimus task level outputs to Multiome pipeline level outputs

# 3.0.5 
2024-01-18 (Date of Last Commit)

* Increased memory for MergeStarOutputs in StarAlign.wdl, RunEmptyDrops in RunEmptyDrops.wdl, OptimusH5ad in H5adUtils.wdl and GeneMetrics in Metrics.wdl
* Added the --soloMultiMappers flag as an optional input to the StarSoloFastq task in the StarAlign.wdl
* Added a check of read2 length to the paired-tag pipeline; this does not affect the Multiome workflow

# 3.0.4
2024-01-05 (Date of Last Commit)

* Added new functionality to the ATAC workflow for paired-tag data, including the option for SnapATAC to pull cell barcodes from the BB tag of the BAM
* Modified the STARsoloFastq task in the StarAlign.wdl so STARsolo can run different types of alignments in a single STARsolo command depending on the counting_mode

# 3.0.3
2023-12-20 (Date of Last Commit)

* Added updated docker to BWAPairedEndAlignment ATAC task to use updated code for distributed bwa-mem2 from Intel
* Removed MergedBAM ATAC and moved BWAPairedEndAlignment ATAC outside of the for loop
* Changed CPU platform to Ice Lake for BWAPairedEndAlignment ATAC task
* Added input parameter input_output_parameter to the Multiome ATAC wdl

# 3.0.2
2023-12-20 (Date of Last Commit)

* JoinMultiomeBarcodes now has dynamic memory and disk allocation

# 3.0.1
2023-12-12 (Date of Last Commit)

* ValidateVcfs now has optional memory parameter; this does not affect this pipeline
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

