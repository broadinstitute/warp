# 3.6.3
2025-06-20 (Date of Last Commit)
* Added reference genome/GTF headers to fragment file via new string inputs; this change does not affect this pipeline

# 3.6.2
2025-06-09 (Date of Last Commit)
* Increased the ulimit in the following tasks: CalculateCellMetrics, CalculateGeneMetrics, CalculateUMIsMetrics; this does not affect the outputs of the pipeline

# 3.6.1
2025-05-27 (Date of Last Commit)
* Increased the ulimit in the STARsoloFastq task in the StarAlign.wdl to 10000; this does not affect the outputs of the pipeline

# 3.6.0
2025-04-02 (Date of Last Commit)
* Removed MergeStarOutput task and updated docker image in alignemnt step (STARsoloFastq) in Optimus; this does not affect the outputs of the pipeline

# 3.5.0
2025-02-25 (Date of Last Commit)
* Refactored the STAR alignment step (STARsoloFastq) in Optimus and removed tasks FastqProcessing and MergeSortBamFiles; we are no longer sharding. We are now running one instance of STAR; this does not affect the outputs of the pipeline

# 3.4.9
2025-02-25 (Date of Last Commit)
* Updated the warp-tools docker image to include an update to the GroupQCs function in sctools; this does not affect the outputs of the pipeline
* Added reference information to the BAM header for Optimus and ATAC workflows; this does not impact Slideseq

# 3.4.8
2025-01-13 (Date of Last Commit)

*  Added a boolean variable is_slidetags; this does not affect the outputs of the pipeline
* Added reference_gtf_file to the output h5ad unstructured metadata

# 3.4.7
2024-12-3 (Date of Last Commit)

* Fixed a bug in the StarSoloFastq task that caused the pipeline to not output a UniqueAndMult-Uniform.mtx when --soloMultiMappers Uniform was passed to STAR; this does not affect the outputs of the pipeline

# 3.4.6
2024-11-15 (Date of Last Commit)

* Added bam validation in the StarSoloFastq task; this does not affect the outputs of the pipeline
* Updated the warp-tools docker; this update changes the way gene_names are identified when creating gene expression h5ad files

# 3.4.5
2024-11-12 (Date of Last Commit)

* Added memory and disk updates to Multiome JoinBarcodes; this does not impact the SlideSeq workflow

# 3.4.4
2024-11-04 (Date of Last Commit)

* Updated the tabix flag in JoinMultiomeBarcodes task in H5adUtils.wdl to use CSI instead of TBI indexing, which supports chromosomes larger than 512 Mbp; this task should not affect the Slide-seq pipeline


# 3.4.3
2024-10-24 (Date of Last Commit)

* Updated the h5adUtils WDL to rename the gene expression library CSV filename; this does not impact slideseq
* Updated the ATAC fragment file output so that it is bgzipped; this does not impact the slideseq workflow
* Updated memory settings for PairedTag; does not impact the Slideseq workflow

# 3.4.2
2024-09-24 (Date of Last Commit)

* Added a python implementation of DoubletFinder to calculate doublet scores in gene expression data; this does not impact the slideseq workflow

# 3.4.1
2024-09-11 (Date of Last Commit)

* Updated warp-tools docker which added create_h5ad_snss2.py to the docker image. This change does not affect the SlideSeq pipeline

# 3.4.0
2024-08-06 (Date of Last Commit)

* Updated the warp-tools docker to calculate mitochondrial reads from unique reads in cell and gene metrics; these metrics are in the cell and gene metrics CSV as well as h5ad


# 3.3.1
2024-08-02 (Date of Last Commit)

* The ubuntu_16_0_4 docker image version was pinned instead of using the latest tag; this does not affect the outputs of the pipeline

# 3.3.0
2024-07-25 (Date of Last Commit)

* Updated the warp-tools docker image to add TSO metrics to the output h5ad and metric CSV files


# 3.2.0
2024-07-11 (Date of Last Commit)

* Updated the Optimus.wdl to run on Azure; cloud_provider is a new, required input

# 3.1.8
2024-07-09 (Date of Last Commit)

* Added new optional input parameter of gex_nhash_id to the STARAlign task; this does not impact the SlideSeq workflow 

# 3.1.7
2024-06-28 (Date of Last Commit)

* Updated the STARsolo parameters for estimating cells to Emptydrops_CR; this does not affect the slideseq pipeline

# 3.1.6
2024-05-20 (Date of Last Commit)

* Updated SnapATAC2 docker to SnapATAC2 v2.6.3; this does not impact the SlideSeq workflow

# 3.1.5
2024-04-12 (Date of Last Commit)

* Updated the input parameters for STARsolo in STARsoloFastq task. These include the parameters: soloCBmatchWLtype, soloUMIdedup and soloUMIfiltering

# 3.1.4
2024-03-26 (Date of Last Commit)

* Updated the median umi per cell metric for STARsolo library-level metrics

# 3.1.3
2024-03-15 (Date of Last Commit)

* Added cell metrics to the library-level metrics CSV; this does not impact the slide-seq pipeline

* Updated the docker for the MergeStarOutput task to include STARsolo v2.7.11a and custom scripts to create a uniform matrix file and scripts to collect library-level metrics from STARsolo output

* Modified the MergeStarOutput to call a custom script for creating a uniform matrix file (mtx) from individual shard mtx files and to create a filtered matrix from the uniform matrix with STARsolo
# 3.1.2

2024-02-28 (Date of Last Commit)

* Updated the Optimus workflow to produce a library-level metrics CSV; this does not impact the slide-seq pipeline

# 3.1.1
2024-02-29 (Date of Last Commit)
* Added mem and disk to inputs of Join Barcodes task of Multiome workflow; does not impact the Slideseq workflow


# 3.1.0
2024-02-07 (Date of Last Commit)

* Updated StarAlign output metrics to include shard ids

# 3.0.1
2024-02-13 (Date of Last Commit)

* Updated the Metrics tasks to exclude mitochondrial genes from reads_mapped_uniquely, reads_mapped_multiple and reads_mapped_exonic, reads_mapped_exonic_as and reads_mapped_intergenic; this does affect the SlideSeq workflow

# 3.0.0
2024-02-12 (Date of Last Commit)

* Updated the SlideSeq WDL output to utilize the h5ad format in place of Loom


# 2.1.6
2024-01-30 (Date of Last Commit)

* Added task GetNumSplits before FastqProcess ATAC task to determine the number of splits based on the bwa-mem2 machine specs; this does affect the SlideSeq workflow
* Added an error message to the BWAPairedEndAlignment ATAC task to ensure that the number of splits equals the number of ranks; this does affect the SlideSeq workflow
* Added an error message to the BWAPairedEndAlignment ATAC task to ensure that the number of R1s equals to the number of R3s; this does affect the SlideSeq workflow

# 2.1.5
2024-01-11 (Date of Last Commit)

* Increased memory for MergeStarOutputs in StarAlign.wdl, RunEmptyDrops in RunEmptyDrops.wdl, OptimusH5ad in H5adUtils.wdl and GeneMetrics in Metrics.wdl
* Added the --soloMultiMappers flag as an optional input to the StarSoloFastq task in the StarAlign.wdl; this does affect the SlideSeq workflow
  
# 2.1.4
2024-01-05 (Date of Last Commit)

* Modified the STARsoloFastq task in the StarAlign.wdl so STARsolo can run different types of alignments in a single STARsolo command depending on the counting_mode; this does affect the SlideSeq workflow

# 2.1.3
2023-12-17 (Date of Last Commit)

* Updated the ATAC WDL for the Multiome BWAPairedEndAlignment and MergedBAM tasks; this does affect the SlideSeq workflow
  

# 2.1.2
2023-11-21 (Date of Last Commit)

* Added the latest warp-tools docker to tasks in the Metrics, FastqProcessing and H5adUtils wdls; this incorporates new input parameter for number of output fastq files to fastqprocess

# 2.1.1
2023-11-20 (Date of Last Commit)

* Added the latest warp-tools docker to the Metrics task; this allows use of REFSEQ references

# 2.1.0
2023-11-03 (Date of Last Commit)

* Updated the Metrics task so that Cell Metrics and Gene Metrics now calculate intronic, intronic_as, exonic, exonic_as, and intergenic metrics from unique reads only using the NH:i:1 tag in the BAM

# 2.0.1
2023-09-21 (Date of Last Commit)
* Added dynamic barcode orientation selection to the ATAC workflow FastqProcess task; this does not impact Slideseq

# 2.0.0
2023-08-22 (Date of Last Commit)

* Updated Slideseq pipeline to include STARsolo v2.7.11a; this impacts the UMI metrics CSV for unassigned genes
* Added sF tag to STARsolo aligner parameters
* Updated TagSort tool for Metrics task to calculate metrics based on the sF tag
* Modified H5adUtils task to include new metrics in the final Optimus h5ad; does not impact slide-seq
* Removed the Dropseq metrics task; this change does not impact the Slideseq workflow

# 1.0.10
2023-07-18 (Date of Last Commit)

* Added STARsolo v2.7.10b metric outputs as an optional pipeline output and an output of the STARalign and MergeSTAR tasks. This does not impact the Slideseq pipeline

# 1.0.9
2023-06-14 (Date of Last Commit)

* Updated the Metrics task WDL for adding a Dropseq task to Optimus; this has no impact on the Slide-seq workflow
* Updated the FastqProcessing.wdl for ATAC. This update has no impact on the SlideSeq workflow
* Updated STARsolo version to v2.7.10b for the StarsoloFastq task; this change does not impact this workflow
* Updated STARsolo argument for counting mode to GeneFull_Ex50pAS; this change does not impact this workflow 


# 1.0.8
2023-05-31 (Date of Last Commit)

* Updated the FastqProcessing.wdl. This update has no impact on the SlideSeq workflow

# 1.0.7
2023-05-11 (Date of Last Commit)

* Updated Docker image for FastqProcessing task to latest warp-tools

# 1.0.6
2023-05-04 (Date of Last Commit)

* Updated the warp-tools docker for Optimus fastqprocess; change does not impact this workflow

# 1.0.5
2023-04-23 (Date of Last Commit)

* Updated the STARalign task; does not affect this workflow

# 1.0.3
2023-04-19

* Updated warp-tools docker which included a fix for a small bug in create_snrna_optimus.py that was causing the script not to run

# 1.0.3
2023-03-27 (Date of Last Commit)
* Removed the following columns from the gene metrics csv and the Loom as the counts were empty/incorrect: reads_unmapped, reads_mapped_exonic, reads_mapped_intronic, reads_mapped_utr, reads_mapped_intergenic, duplicate_reads. We also removed duplicate_reads from the cell metrics csv and the Loom

# 1.0.2
2023-03-15 (Date of Last Commit)

* Removed the following columns from the cell metrics csv and the Loom as the counts were empty/incorrect: reads_unmapped, reads_mapped_exonic, reads_mapped_intronic, reads_mapped_utr, reads_mapped_intergenic
* Updated warp-tools docker. The latest loom building script updates the optimus_output_schema_version from 1.0.0 to 1.0.1 to capture the metrics changes listed above 


# 1.0.1
2023-02-28 (Date of Last Commit)

* Added a new task to the workflow that reads the tar_star_reference file to obtain the genomic reference source, build version, and annotation version and outputs the information as a txt file

# 1.0.0
2023-02-13 (Date of Last Commit)

* The first major version release for the SlideSeq pipeline
