# 2.2.3
2026-01-21 (Date of Last Commit)

* Moved inputs into new Google buckets. This change does not affect the outputs of the pipeline

# 2.2.2
 2025-06-20 (Date of Last Commit)

* Added reference genome/GTF headers to fragment file via new string inputs; this change does not affect this pipeline

# 2.2.1
2025-05-27 (Date of Last Commit)

* Increased the ulimit in the STARsoloFastq task in the StarAlign.wdl to 10000; this does not affect the outputs of the  pipeline

# 2.2.0
2025-04-02 (Date of Last Commit)

* Removed MergeStarOutput task and updated docker image in alignemnt step (STARsoloFastq) in Optimus; this does not affect the outputs of the pipeline

# 2.1.0
2025-03-19 (Date of Last Commit)

* Refactored the STAR alignment step (STARsoloFastq) in Optimus and removed tasks FastqProcessing and MergeSortBamFiles; we are no longer sharding. We are now running one instance of STAR; this does not affect the outputs of the pipeline

# 2.0.8
2025-02-25 (Date of Last Commit)

* Updated the warp-tools docker image to include an update to the GroupQCs function in sctools; this does not affect the outputs of the pipeline
* Added reference information to the BAM header for Optimus and ATAC; does not impact snSS2

# 2.0.7
2025-01-13 (Date of Last Commit)

*  Added a boolean variable is_slidetags; this does not affect the outputs of the pipeline
* Added reference_gtf_file to the output h5ad unstructured metadata

# 2.0.6
2024-11-15 (Date of Last Commit)

* Fixed a bug in the StarSoloFastq task that caused the pipeline to not output a UniqueAndMult-Uniform.mtx when --soloMultiMappers Uniform was passed to STAR; this does not affect the outputs of the pipeline

# 2.0.5
2024-11-15 (Date of Last Commit)

* Added bam validation in the StarSoloFastq task; this does not affect the outputs of the pipeline

# 2.0.4
2024-11-12 (Date of Last Commit)

* Added memory and disk updates to Multiome JoinBarcodes; this does not impact the snSS2 workflow

# 2.0.3
2024-11-04 (Date of Last Commit)

* Updated the tabix flag in JoinMultiomeBarcodes task in H5adUtils.wdl to use CSI instead of TBI indexing, which supports chromosomes larger than 512 Mbp; this task should not affect the snSS2 pipeline


# 2.0.2
2024-10-23 (Date of Last Commit)

* Updated the h5adUtils WDL to rename the gene expression library CSV filename; this does not impact slideseq
* Updated the ATAC fragment file output so that it is bgzipped; this does not impact the Multi-snSS2 workflow
* Updated memory settings for PairedTag; does not impact the snSS2 workflow

# 2.0.1
2024-09-24 (Date of Last Commit)

* Added a python implementation of DoubletFinder to calculate doublet scores in gene expression data; this does not affect the snSS2 workflow

# 2.0.0
2024-09-11 (Dat of Last Commit)

* Added h5ad as a format option for the cell by gene matrix output. The h5ad has the same layers and global attributes (unstructured data in h5ad) as the previous Loom output

# 1.4.2
2024-08-25-02 (Dat of Last Commit)

* The ubuntu_16_0_4 docker image version was pinned instead of using the latest tag; this does not affect the outputs of the pipeline

# 1.4.1
2024-07-25 (Dat of Last Commit)

* Updated the warp-tools docker image to add TSO metrics to the output h5ad and metric CSV files; this does not impact the snSS2 workflow

# 1.4.0
2024-07-11 (Date of Last Commit)

* Updated the PairedTag.wdl to run on Azure. cloud_provider is a new, required input.
* Added new optional input parameter of gex_nhash_id to the STARAlign task; this does not impact the MultiSampleSmartSeq2SingleNucleus workflow 

# 1.3.5
2024-06-28 (Date of Last Commit)

* Updated the STARsolo parameters for estimating cells to Emptydrops_CR; this does not impact the snSS2 pipeline

# 1.3.4
2024-04-12 (Date of Last Commit)

* Updated the input parameters for STARsolo in STARsoloFastq task. These include the parameters: soloCBmatchWLtype, soloUMIdedup and soloUMIfiltering

# 1.3.3
2024-03-26 (Date of Last Commit)

* Updated the median umi per cell metric for STARsolo library-level metrics

# 1.3.2
2024-03-15 (Date of Last Commit)

* Added cell metrics to the library-level metrics CSV; this does not impact the Single-nucleus Multi Sample Smartseq pipeline

* Updated the docker for the MergeStarOutput task to include STARsolo v2.7.11a and custom scripts to create a uniform matrix file and scripts to collect library-level metrics from STARsolo output

* Modified the MergeStarOutput to call a custom script for creating a uniform matrix file (mtx) from individual shard mtx files and to create a filtered matrix from the uniform matrix with STARsolo
# 1.3.1
2024-02-28 (Date of Last Commit)

* Updated the Optimus workflow to produce a library-level metrics CSV; this does not impact the Single-nucleus Multi Sample Smart-seq2 pipeline

# 1.3.0
2024-01-22 (Date of Last Commit)

* Updated StarAlign output metrics to include shard ids
 
 # 1.2.28
2024-01-11 (Date of Last Commit)

* Increased memory for MergeStarOutputs in StarAlign.wdl, RunEmptyDrops in RunEmptyDrops.wdl, OptimusH5ad in H5adUtils.wdl and GeneMetrics in Metrics.wdl
* Added the --soloMultiMappers flag as an optional input to the StarSoloFastq task in the StarAlign.wdl; this does affect the MultiSampleSmartSeq2SingleNucleus workflow
 
# 1.2.27
2024-01-05 (Date of Last Commit)

* Modified the STARsoloFastq task in the StarAlign.wdl so STARsolo can run different types of alignments in a single STARsolo command depending on the counting_mode; this does affect the MultiSampleSmartSeq2SingleNucleus workflow

# 1.2.26
2023-08-22 (Date of Last Commit)

* Updated Optimus pipeline to include STARsolo v2.7.11a; does not impact snSS2
* Added sF tag to STARsolo aligner parameters; does not impact snSS2
* Updated TagSort tool for Optimus Metrics task to calculate metrics based on the sF tag; does not impact snSS2
* Modified H5adUtils task to include new metrics in the final Optimus h5ad; does not impact snSS2

# 1.2.25
2023-07-18 (Date of Last Commit)

* Added STARsolo v2.7.10b metric outputs as an optional pipeline output and an output of the STARalign and MergeSTAR tasks. This does not impact the snSS2 pipeline
* Updated the CountAlignments task in the FeatureCounts.wdl to use a new docker image. This change does not affect the MultiSampleSmartSeq2SingleNucleus pipeline


# 1.2.24
2023-06-23 (Date of Last Commit)

* Added a move command to work around a maximum file path length in featureCounts
* Updated STARsolo version to v2.7.10b for the StarsoloFastq task; does not impact this workflow
* Updated STARsolo argument for counting mode to GeneFull_Ex50pAS; does not impact this workflow

# 1.2.23 
2023-05-04 (Date of Last Commit)

* Updated the CheckInputs WDL for the Optimus workflow. This changes does impact snSS2

# 1.2.22

2023-04-23 (Date of Last Commit)

* Updated the STARalign task; does not affect this workflow

# 1.2.21
2023-04-19 (Date of Last Commit)

* Updated warp-tools docker which included a fix for a small bug in create_snrna_optimus.py that was causing the script not to run


# 1.2.20
2023-03-27 (Date of Last Commit)

* SlideSeq-specific and Optimus-specific changes to Metrics.wdl. This change does not affect the MultiSampleSmartSeq2SingleNucleus pipeline

# 1.2.19
2023-03-15 (Date of Last Commit)

* SlideSeq-specific and Optimus-specific changes to Metrics.wdl. This change does not affect the MultiSampleSmartSeq2SingleNucleus pipeline
* Updated warp-tools docker to support the Optimus and SlideSeq changes 

# 1.2.18
2023-02-28 (Date of Last Commit)

* Added a new task to the workflow that reads the tar_star_reference file to obtain the genomic reference source, build version, and annotation version and outputs the information as a txt file

# 1.2.17
2023-02-13 (Date of Last Commit)

* SlideSeq-specific changes to LoomUtils.wdl. This change does not affect the MultiSampleSmartSeq2SingleNucleus pipeline

# 1.2.16
2023-02-07 (Date of Last Commit)

* Updated the input checks for the Optimus pipeline task; this has no effect on this pipeline
* Added disk to the checkInputArrays task 

# 1.2.15
2023-01-23 (Date of Last Commit)

* Added 'Disk' to task runtime sections to support running on Azure
* Addressed mb/gb memory specification inconsistencies in LoomUtils and CheckInput

# 1.2.14
2022-09-20 (Date of Last Commit)

* Updated tasks in StarAlign.wdl to use an updated STAR docker image.

# 1.2.13

2022-09-01 (Date of Last Commit)

* Updated CheckInputs.wdl to use a lightweight alpine-bash image.

# 1.2.12
2022-08-31 (Date of Last Commit)

* Updated CountAlignments to use an updated docker image.

# 1.2.11
2022-08-23 (Date of Last Commit)

* Removed an unused script in pytools docker image.

# 1.2.10
2022-08-16 (Date of Last Commit)

* Updated LoomUtils.wdl to use a consolidated python utilities docker image. This change does not affect the MultiSampleSmartSeq2SingleNucleus pipeline.

# 1.2.9
2022-08-08 (Date of Last Commit)

* Updated TrimAdapters runtime docker URL.

# 1.2.8
2022-07-21 (Date of Last Commit)

* Updated STARsoloFastq runtime docker URL.

# 1.2.7
2022-05-18 (Date of Last Commit)

* Updated merge npz docker in StarAlign.wdl

# 1.2.6
2022-04-22 (Date of Last Commit)

* Updated LoomUtils.wdl for a task in the Optimus pipeline. This change does not affect the MultiSampleSmartSeq2SingleNucleus pipeline.

# 1.2.5
2022-04-19 (Date of Last Commit)

* Updated to Picard version 2.26.10 and GATK version 4.2.6.1 to address log4j vulnerabilities

# 1.2.4

2022-04-06 (Date of Last Commit)

* Updated STARsoloFastq task in StarAlign.wdl; this change does not affect the Mutl-snSS2 workflow.

# 1.2.3
2022-02-25 (Date of Last Commit)

* Updated LoomUtils.wdl for a task in the Optimus pipeline. This change does not affect the MultiSampleSmartSeq2SingleNucleus pipeline.

# 1.2.2
2022-02-10 (Date of Last Commit)

* Rebuilt a docker to merge outputs of STAR in StarAlign.wdl task and moved it to a public location.

# 1.2.1
2022-02-07 (Date of Last Commit)

* Updated a task in STARalign.wdl related to the Optimus pipeline. This pipeline has not been changed.

# 1.2.0
2022-01-07 (Date of Last Commit)

* Fixed missing metadata issue in the loom file
# 1.1.2
2021-11-19 (Date of Last Commit)

* Updated STARsoloFastq to use 'output_bam_basename' to name the aligned bam. This is consistent with versions 4.2.7 and older. This change has no impact MultiSampleSmartSeq2SingleNucleus
# 1.1.1
2021-11-15 (Date of Last Commit)
* Updated remove-reads-on-junctions.py in the FeatureCounts.wdl to use python3 instead of python2.

# 1.1.0
2021-09-16 (Date of Last Commit)

* Removed the Smart-seq2 Single Nucleus workflow (SmartSeq2SingleNucleus.wdl) and changed the workflow tasks to run multiple samples in the same VM. This change is expected to make the pipeline faster and cheaper.
* Renamed the StarAlignFastq.StarAlignFastqPairedEnd task to StarAlign.StarAlignFastqMultisample

# 1.0.4
2021-09-10 (Date of Last Commit)

Added the option "--soloBarcodeReadLength 0" STARsoloFastq task to support alignment in Optimus. This change has no impact on MultiSampleSmartSeq2SingleNucleus

# 1.0.3
2021-09-02 (Date of Last Commit)

* Added a new StarSolo task for Optimus in the StarAlign.wdl. However, the same wdl
  contains other Star tasks that are used in the smartseq2 single nuclei for paired and 
  single stranged fastq files. As a result, the smartseq2 processing is not expected to 
  change. 

# 1.0.2
2021-08-02 (Date of Last Commit)

* Increased the version number to make new release tag for Dockstore 

# 1.0.1
2021-07-19 (Date of Last Commit)

* Updated SmartSeq2 to accommodate spaces in input_name

# 1.0.0

2021-05-17 (Date of First Commit)

* This is the first release of the Smart-seq2 Multi-Sample Single Nuclei workflow
