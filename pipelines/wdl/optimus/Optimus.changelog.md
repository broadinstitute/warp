# 8.0.4
2026-01-22 (Date of Last Commit)

* Added a new, defaulted input cellbender_memory_GB to Optimus; this does not affect the outputs of the pipeline

# 8.0.4
25-07-31 (Date of Last Commit)

* Reorganized all WDL pipelines into the wdl pipeline directory

# 8.0.3
2025-06-20 (Date of Last Commit)

* Added reference genome/GTF headers to fragment file via new string inputs; this change does not affect this pipeline

# 8.0.2
2025-06-09 (Date of Last Commit)

* Removed quotes from bootDiskSizeGb in RunEmptyDrops task to be compatible with Google Batch; this does not affect the outputs of the pipeline
* Increased the ulimit in the following tasks: CalculateCellMetrics, CalculateGeneMetrics, CalculateUMIsMetrics; this does not affect the outputs of the pipeline

# 8.0.1
2025-05-27 (Date of Last Commit)

* Increased the ulimit in the STARsoloFastq task in the StarAlign.wdl to 10000; this does not affect the outputs of the pipeline

# 8.0.0
2025-04-02 (Date of Last Commit)

* Implemented a unified STARsolo execution strategy to ensure consistent and accurate cell barcode correction across the entire dataset. This update resolves discrepancies that previously arose from sharded (partitioned) processing, where each shard independently corrected barcodes using incomplete local priors. By consolidating barcode frequency calculations and applying correction globally, the pipeline now mirrors the behavior of DropSeq and Cell Ranger
* Removed boolean variable is_slidetags; no longer needed with new updates
* Refactored the STAR alignment step and removed tasks FastqProcessing and MergeSortBamFiles
* Added parameters for STARsoloFastq task, including cpu_platform_star, mem_size_star, cpu_star, disk_star, limitBAMsortRAM_star, and outBAMsortingBinsN_star, for dynamic allocation of resources depending on input size
* Removed MergeStarOutput tasks; added necessary parts of MergeStarOutput task to the STAR alignment step (STARsoloFastq). Additional outputs added to STARsoloFastq task as a result; this includes row_index, col_index, sparse_counts, library_metrics, mtx_files, filtered_mtx_files and cell_reads_out
* Updated the STAR docker image to include Samtools and Python

# 7.9.2
2025-02-25 (Date of Last Commit)

* Updated the warp-tools docker image to include an update to the GroupQCs function in sctools; this does not affect the outputs of the pipeline
* Added reference information to the BAM header

# 7.9.1
2025-01-13 (Date of Last Commit)

* Added a boolean variable is_slidetags; set to false by default, but set to true if the Slide-Tags pipeline is calling Optimus
* Added reference_gtf_file to the output h5ad unstructured metadata

# 7.9.0
2024-12-05 (Date of Last Commit)

* Added an optional task to the Optimus.wdl that will run CellBender on the Optimus output h5ad file

# 7.8.4
2024-12-3 (Date of Last Commit)

* Fixed a bug in the StarSoloFastq task that caused the pipeline to not output a UniqueAndMult-Uniform.mtx when --soloMultiMappers Uniform was passed to STAR

# 7.8.3
2024-11-22 (Date of Last Commit)

* Added bam validation in the StarSoloFastq task; this does not affect the outputs of the pipeline
* Updated the warp-tools docker; this update changes the way gene_names are identified when creating gene expression h5ad files

# 7.8.2
2024-11-12 (Date of Last Commit)

* Added memory and disk updates to Multiome JoinBarcodes; this does not impact the Optimus workflow

# 7.8.1
2024-11-04 (Date of Last Commit)

* Updated the tabix flag in JoinMultiomeBarcodes task in H5adUtils.wdl to use CSI instead of TBI indexing, which supports chromosomes larger than 512 Mbp; this task should not affect the Optimus pipeline


# 7.8.0
2024-10-23 (Date of Last Commit)

* Renamed the input expected_cells to gex_expected_cells
* Updated gex_expected_cells to a required output
* Reformatted the library CSV output filename to remove an extra gex
* Updated the ATAC fragment file output so that it is bgzipped; this does  not impact the Optimus workflow
* Updated memory settings for PairedTag; does not impact the Optimus workflow

# 7.7.0
2024-09-24 (Date of Last Commit)

* Added a python implementation of DoubletFinder to calculate doublet scores in gene expression data; percent doublets are now available as a library-level metric and individual doublet scores for cell barcodes are in the h5ad
* Updated gene_names in the final h5ad to be unique

# 7.6.1
2024-09-11 (Date of Last Commit)

* Updated warp-tools docker which added create_h5ad_snss2.py to the docker image. This change does not affect the Optimus pipeline

# 7.6.0
2024-08-06 (Date of Last Commit)

* Updated the warp-tools docker to calculate mitochondrial reads from unique reads in cell and gene metrics; these metrics are in the cell and gene metrics CSV as well as h5ad

# 7.5.1
2024-08-02 (Date of Last Commit)

* The ubuntu_16_0_4 docker image version was pinned instead of using the latest tag; this does not affect the outputs of the pipeline

# 7.5.0
2024-07-25 (Date of Last Commit)

* Updated the warp-tools docker image to add TSO metrics to the output h5ad and metric CSV files
* Update the library-level metrics to include new TSO metrics and NHashID descriptor

# 7.4.0
2024-07-11 (Date of Last Commit)

* Updated the Optimus.wdl to run on Azure. cloud_provider is a new, required input.
* Updated GermlineVariantDiscovery, BamProcessing, DragenTasks, Qc, and Utilities tasks to allow multi-cloud dockers

# 7.3.0
2024-07-09 (Date of Last Commit)

* Added new optional input parameter of gex_nhash_id, a string identifier for a library aliquot that is echoed in the h5ad cell by gene matrix (in the data.uns) and the library metrics CSV output; default is set to null 

# 7.2.0
2024-06-28 (Date of Last Commit)

* Updated the STARsolo parameters for estimating cells to Emptydrops_CR
* Added an optional input for expected cells which is used for metric calculation

# 7.1.0
2024-05-20 (Date of Last Commit)

* Updated SnapATAC2 docker to SnapATAC2 v2.6.3; this does not impact the Optimus workflow

# 7.0.0
2024-04-24 (Date of Last Commit)

* Updated the input parameters for STARsolo in STARsoloFastq task. These include the parameters: soloCBmatchWLtype, soloUMIdedup and soloUMIfiltering
* Added "Uniform" as the default string for STARsolo multimapping parameters

# 6.6.1
2024-03-26 (Date of Last Commit)

* Updated the median umi per cell metric for STARsolo library-level metrics

# 6.6.0
2024-03-15 (Date of Last Commit)

* Added cell metrics to the library-level metrics

* Updated the docker for the MergeStarOutput task to include STARsolo v2.7.11a and custom scripts to create a uniform matrix file and scripts to collect library-level metrics from STARsolo output

* Modified the MergeStarOutput to call a custom script for creating a uniform matrix file (mtx) from individual shard mtx files and to create a filtered matrix from the uniform matrix with STARsolo

# 6.5.0
2024-02-28 (Date of Last Commit)

* Added a library-level metrics CSV as output of the Optimus workflow; this iteration includes read-level metrics 

# 6.4.1
2024-02-29 (Date of Last Commit)
* Added mem and disk to inputs of Join Barcodes task of Multiome workflow; does not impact the Optimus workflow


# 6.4.0
2024-02-21 (Date of Last Commit)
* Updated StarAlign.MergeStarOutput to add a shard number to the metrics files
* Removed ref_genome_fasta input from Optimus WDL and JSON

# 6.3.6
2024-02-07 (Date of Last Commit)
* Updated the Metrics tasks to exclude mitochondrial genes from reads_mapped_uniquely, reads_mapped_multiple and reads_mapped_exonic, reads_mapped_exonic_as and reads_mapped_intergenic

# 6.3.5
2024-01-30 (Date of Last Commit)
* Added task GetNumSplits before FastqProcess ATAC task to determine the number of splits based on the bwa-mem2 machine specs
* Added an error message to the BWAPairedEndAlignment ATAC task to ensure that the number of splits equals the number of ranks
* Added an error message to the BWAPairedEndAlignment ATAC task to ensure that the number of R1s equals to the number of R3s

# 6.3.4
2024-01-11 (Date of Last Commit)
* Increased memory for MergeStarOutputs in StarAlign.wdl, RunEmptyDrops in RunEmptyDrops.wdl, OptimusH5ad in H5adUtils.wdl and GeneMetrics in Metrics.wdl
* Added the --soloMultiMappers flag as an optional input to the StarSoloFastq task in the StarAlign.wdl

# 6.3.3
2024-01-05 (Date of Last Commit)
* Modified the STARsoloFastq task in the StarAlign.wdl so STARsolo can run different types of alignments in a single STARsolo command depending on the counting_mode

# 6.3.2
2023-12-20 (Date of Last Commit)
* Updated the ATAC WDL for the Multiome BWAPairedEndAlignment and MergedBAM tasks; this does affect the Optimus workflow
  
# 6.3.1
2023-12-20 (Date of Last Commit)
* JoinMultiomeBarcodes now has dynamic memory and disk allocation; this does affect the Optimus workflow

# 6.3.0
2023-12-04 (Date of Last Commit)

* Updated the h5ad utils WDL for the Multiome JoinBarcodes task; this does affect the Optimus workflow

# 6.2.2
2023-11-29 (Date of Last Commit)

* Added the latest warp-tools docker to tasks in the Metrics, FastqProcessing and H5adUtils wdls; this incorporates new input parameter for number of output fastq files to fastqprocess and allows use of REFSEQ references

# 6.2.0
2023-11-03 (Date of Last Commit)

* Updated the Metrics task so that Cell Metrics and Gene Metrics now calculate intronic, intronic_as, exonic, exonic_as, and intergenic metrics from unique reads only using the NH:i:1 tag in the BAM

# 6.1.2
2023-10-20 (Date of Last Commit)

* Removed the dropna from the H5adUtils WDL for the JoinBarcodes task; this change does not impact Optimus outputs
* Added a new task to the H5adUtils WDL to combine Multiome barcodes in h5ad outputs; this does not impact the individual Optimus workflow

# 6.1.0
2023-09-21 (Date of Last Commit)

* Added dynamic barcode orientation selection to ATAC workflow FastqProcess task; this has no impact on Optimus

# 6.0.0
2023-08-22 (Date of Last Commit)

* Updated Optimus pipeline to include STARsolo v2.7.11a
* Added sF tag to STARsolo aligner parameters
* Updated TagSort tool for Optimus Metrics task to calculate metrics based on the sF tag
* Modified H5adUtils task to include new metrics in the final Optimus h5ad
* Removed the Dropseq metrics task

# 5.8.4
2023-07-18 (Date of Last Commit)

* Added STARsolo v2.7.10b metric outputs as an optional pipeline output and an output of the STARalign and MergeSTAR tasks

# 5.8.3
2023-06-23 (Date of Last Commit)

* Updated STARsolo version to v2.7.10b for the StarsoloFastq task
* Updated STARsolo argument for counting mode to GeneFull_Ex50pAS 
* Updated the FastqProcessing.wdl. This update has no impact on the Optimus workflow
* Added h5ad as a format option for the cell by gene matrix output. The h5ad has the same layers and global attributes (unstructured data in h5ad) as the previous Loom output
* Added Dropseq cell metrics to Multiome and Optimus workflows
* Updated the FastqProcessing.wdl for ATAC. This update has no impact on the Optimus workflow 


# 5.8.2
2023-05-11 (Date of Last Commit)

* Updated the Docker image for Slideseq FastqProcessing. This update has no impact on the Optimus workflow

# 5.8.1
2023-05-04 (Date of Last Commit)

* Updated inputs for the FastqProcessing task, which now requires read structure. This is dynamically calculated from the CheckInputs task

# 5.8.0
2023-04-24 (Date of Last Commit)

* Modified the stranded input parameter to be called star_strand_mode; the default is now set to Forward and the other options include Unstranded and Reverse

# 5.7.5
2023-04-19 (Date of Last Commit)

* Updated warp-tools docker which included a fix for a small bug in create_snrna_optimus.py that was causing the script not to run


# 5.7.4
2023-03-27 (Date of Last Commit)

* Removed the following columns from the gene metrics csv and the Loom as the counts were empty/incorrect: reads_unmapped, reads_mapped_exonic, reads_mapped_intronic, reads_mapped_utr, reads_mapped_intergenic, duplicate_reads. We also removed duplicate_reads from the cell metrics csv and the Loom


# 5.7.3
2023-03-15 (Date of Last Commit)

* Removed the following columns from the cell metrics csv and the Loom as the counts were empty/incorrect: reads_unmapped, reads_mapped_exonic, reads_mapped_intronic, reads_mapped_utr, reads_mapped_intergenic
* Updated warp-tools docker. The latest loom building script updates the optimus_output_schema_version from 1.0.0 to 1.0.1 to capture the metrics changes listed above 

# 5.7.2
2023-02-28 (Date of Last Commit)

* Added a new task to the worklow that reads the tar_star_reference file to obtain the genomic reference source, build version, and annotation version and outputs the information as txt file.

# 5.7.1
2023-02-13 (Date of Last Commit)

* SlideSeq-specific changes to FastqProcessing.wdl, LoomUtils.wdl, Metrics.wdl, and StarAlign.wdl. This change does not affect the Optimus pipeline.

# 5.7.0
2023-02-16 (Date of Last Commit)

* Changed the chemistry input to an int tenx_chemistry_version that accepts either 2 or 3.
* Added new whitelist_v2 and whitelist_v3 inputs that are set to a public references and selected based on chemistry.
* Added new checks to the checkOptimusInput task to verify that the chemistry is either v2 or v3.
* Added new checks to the checkOptimusInput that verify the read1 FASTQ is the correct chemistry based on read length; these checks can be ignored with new boolean input ignore_r1_read_length. 
* Updated warp-tools docker version in FastqProcessing.wdl and Metrics.wdl due to previous bug fix (broadinstitute/warp-tools#18). The issue was integer division that caused the percent of mitochondrial molecules to always be calculated as zero.
* Dynamically sized disk in the checkOptimusInput task 

# 5.6.2
2023-02-07 (Date of Last Commit)

* Prepended the input_id to the name of the output file in both the CalculateCellMetrics and CalculateGeneMetrics tasks in the Metrics.wdl.

# 5.6.1
2023-01-23 (Date of Last Commit)

* Added 'Disk' to task runtime sections to support running on Azure
* Updated the emptyDrops container to address concerns outlined in #772 - avoiding the usage of the root directory inside the container. This also includes some optimizations: moved image to GCR instead of Quay, conformed to (most) of our docker style guideline, build time decreased to 1 hour from 1.5-2 hours, and the image size reduced to 1.5GB from 3GB.
* EmptyDrops container has been upgraded to use R 4.2.2 and BiocManager 3.16
* Addressed mb/gb memory specification inconsistencies in LoomUtils and CheckInput

# 5.6.0
2022-12-06 (Date of Last Commit)

* Updated Metrics.wdl and Optimus.wdl to take an optional inputs for mitochondrial gene names.
* Updated FastqProcessing.wdl and Metrics.wdl to use the warp-tools container.

# 5.5.5
2022-09-20 (Date of Last Commit)

* Updated tasks in StarAlign.wdl to use an updated STAR docker image. 

# 5.5.4
2022-09-01 (Date of Last Commit)

* Updated CheckInputs.wdl to use a lightweight alpine-bash image.

# 5.5.3
2022-08-23 (Date of Last Commit)

* Removed an unused script in pytools docker image and removed unused ConvertStarOutputs task.

# 5.5.2
2022-08-16 (Date of Last Commit)

* Updated LoomUtils.wdl and StarAlign.wdl to use a rebuilt python utilities docker.

# 5.5.1
2022-07-21 (Date of Last Commit)

* Updated STARsoloFastq runtime docker URL.

# 5.5.0
2022-05-18 (Date of Last Commit)

* Updated merge npz docker in StarAlign.wdl to fix a bug in the output loom matrix where gene names were inapporpriately assigned to counts. Any data previously processed with Optimus version 5.0.0 and above should be re-analyzed.
 

# 5.4.3
2022-04-22 (Date of Last Commit)

* Updated Optimus to not run emptydrop step in sn_rna mode.

# 5.4.2
2022-04-21 (Date of Last Commit)

* Updated to Picard version 2.26.10 and GATK version 4.2.6.1 to address log4j vulnerabilities

# 5.4.1
2022-04-21 (Date of Last Commit)

* Fixing syntax in changelog documentation

# 5.4.0
2022-04-06 (Date of Last Commit)

* Updated the STARsoloFastq task in StarAlign.wdl to run STARsolo independently with \"Gene\" COUNTING_MODE when the Optimus input parameter `count_exons` is set to true.
* Updated the MergeStarOutput task in StarAlign.wdl to run an updated script for merging. The previous version included a bug where certain barcodes were getting zero counts after merging.
* Changed the npz output names to include input_id as a prefix. All outputs are prefixed with input_id and are differentiated from run to run with different inputs.

# 5.3.1
2022-04-04 (Date of Last Commit)

* Update task CopyFilesFromCloudToCloud in Utilities.wdl, this update has no effect on this pipeline

# 5.3.0
2022-02-22 (Date of Last Commit)

* Added an optional flag count_exons as the Optimus workflow input with default value of false. 
  If this flag is true, the pipeline adds two layers to the loom file: one for reads aligned to the
  entire gene region, and the second layer is will be a count matrix of reads aligned to only exons.

# 5.2.1
2022-02-10 (Date of Last Commit)

* Rebuilt a docker to merge outputs of STAR in in StarAlign.wdl task and moved it to a public location.

# 5.2.0
2022-01-07 (Date of Last Commit)

* Updated the pipeline to split the fastq files and run parallel STARsolo jobs.
* Added SplitFastq.wdl task to split fastq files by cell barcodes such that each shard gets all reads from the same cell.
* Added MergeSortBam.wdl task to merge bam files from different shards

# 5.1.3
2022-01-07 (Date of Last Commit)

* Updated LoomUtils.wdl to fix a missing metadata issue for the Smart-seq2 Single Nucleus Multi-Sample pipeline. This task update does not affect the Optimus pipeline

# 5.1.2
2021-11-19 (Date of Last Commit)

* Updated STARsoloFastq to use 'output_bam_basename' to name the aligned bam. This is consistent with versions 4.2.7 and older

# 5.1.1
2021-09-13 (Date of Last Commit)

* Updated Picard.wdl and LoomUtils.wdl for Single Nucleus SmartSeq2. These changes do not affect Optimus

# 5.1.0
2021-09-10 (Date of Last Commit)

* Added the option "--soloBarcodeReadLength 0" to STARsoloFastq task to ignore Barcode + UBI read of incorrect length

# 5.0.0
2021-08-30 (Date of Last Commit)

* Replaced STAR alignment with STARsolo and modified the structure of the WDL to utilize the UMI and barcode correction from STARsolo. In the new implementation of Optimus, STARsolo uses the FASTQ file as input and directly creates a count matrix file and a BAM file. No updates have been made to the inputs or outputs. The outputs for this version are identical to the outputs for the previous Optimus version. 
* Updated GoTC base image to AppSec approved 
* Updated BWA version for GoTC image from 0.7.15.r1140 to 0.7.15

# 4.2.7
2021-08-02 (Date of Last Commit)

* Increased the version number to make new release tag for Dockstore 

# 4.2.6
2021-07-19 (Date of Last Commit)

* Updated SmartSeq2 to accommodate spaces in input_name

# 4.2.5

2021-05-24 (Date of Last Commit)

* Updated Picard task to support single nucleus SS2. Changes do not affect Optimus
* We also updated STAR to 2.7.9a

# 4.2.4

2021-04-07 (Date of Last Commit)

* Changed the name of the wdl StarAlignSingleEnd.wdl to StarAlign.wdl
* Added a star dockerfile to STAR version 2.7.8a

# 4.2.3

2021-02-23 (Date of Last Commit)

* Made changes to emptydrops tool wrappper script to not fail in cases with small number of cells, instead, create empty drop result files with NAs.
* Updated the docker in RunEmptyDrops.wdl task to 0.1.4 Updated emptyDropsWrapper.R in the docker

# 4.2.2

2021-01-04 (Date of Last Commit)

* Added an optional input for the pipeline to read in stranded mode which has a default of false

# 4.2.1

2020-12-04 (Date of Last Commit)

* Updated the docker in LoomUtils.wdl task to 0.0.6

# 4.2.0

2020-12-04 (Date of Last Commit)

* Added "Gene" as row attribute in the loom file duplicating "gene_names" to make the output loom compatible with scanpy 
* Added "CellID" (duplicate of ""cell_names") and "input_id" as column attributes in the loom file to make the output loom compatible with scanpy and Cumulus
* Updated the docker in LoomUtils.wdl task to 0.0.5 to incorporate the above changes

# 4.1.8

2020-11-24 (Date of Last Commit)

* Made CPU, memory, and disk optional parameters for all tasks

# 4.1.7

2020-11-05 (Date of Last Commit)

* Increased memory and disk on several tasks

# 4.1.6

2020-11-03 (Date of Last Commit)

* Updated the docker for FastqProcessing task to version v0.3.12. This version of FastqProcessing solves 32 bit unsigned integer overflow error for large files and also disables checking of EOF magic string in BGZF files to accomodate fastq.gz file without this magic string 

# 4.1.5

2020-10-26 (Date of Last Commit)

* Updated the docker in LoomUtils.wdl task to 0.0.4-ss2-loom-fix-1

# 4.1.4

2020-10-22

* Added a check for file naming to FastqProcessing task

# 4.1.3

2020-10-19 (Date of Last Commit)

* Updated memory on fastqprocessing and CalculateGeneMetrics tasks
* Reduced --bam-size to 1 GB on fastqprocessing task

# 4.1.2

2020-10-13 (Date of Last Commit)

* Added a new docker in LoomUtils.wdl

# 4.1.1

2020-10-07 (Date of Last Commit)

* Removed extra trailing slash in ouput directory from cloud to cloud copy job

* Removed fastq_suffix optional input - the pipeline now dynamically determines if a file is zipped

# 4.1.0

2020-10-05 (Date of Last Commit)

* Updated sctools dockers and made them consistent across the Optimus pipeline

# 4.0.2

2020-09-30 (Date of Last Commit)

* Corrected the path to the FastqProcessing WDL


# 4.0.1

2020-09-28 (Date of Last Commit)

* Refactored the pipeline to preprocess fastqs using the task `FastqProcessing`. Outputs are identical and the pipeline should be significantly faster

# 4.0.0

2020-08-10 (Date of Last Commit)
### Breaking changes
* Changed sample_id to input_id
### Non-breaking changes 
* Added input_name as an optional input for user provided sample_id
* Passed pipeline_version to output loom file  
* Added input_id_metadata_field and input_name_metadata_field as optional input


# 3.0.1

2020-07-21 (Date of Last Commit)

* Changed the imports to relative imports to support Dockstore->Terra release
* This version and all future versions have been scientifically validated on 10X 3' snRNAseq data for all supported references

# 3.0.0

2020-06-10 (Date of Last Commit)

* Removed the Zarr formatted matrix and metrics outputs and replaced with Loom
* Removed EmptyDrops for sn_rna mode
* Updated the Loom file attribute names: CellID to cell_names, Gene to gene_names, and Accession to ensembl_ids
* Added metrics for mitochondrial reads
* Added an optional input for the BAM basename; this input is listed as ‘bam_output_basename’and the default is 'sample_id'
* Added a new counting_mode parameter to Optimus workflow which enables processing of single-nuclei datasets
* Updated Drop-seq tools to v2.3.0; this update is only used when the workflow is set to the single-nuclei mode (counting_mode = sn_rna) 
* Updated sctools to support the single-nuclei parameter (counting_mode = sn_rna) 
* Added tests for running the workflow when counting_mode = sn_rna 
* Updated the Loom output to include a global attribute describing the counting mode
* Added new example datasets that can be used with the Optimus workflow
* Updated the README documentation to detail the new counting_mode parameter, describe example datasets, and to include a new FAQ section

# 2.0.0
2020-02-08 (Date of Last Commit)

* Fixed a bug that resulted in emptyDrops output being incorrect
* Updated the workflow to WDL 1.0

# 1.4.0

2019-11-08 (Date of Last Commit)

* Addition of support for V3 chemistry
* Addition of input parameter validation step
* Greatly improved documentation
* Improvements to ZARR output

# 1.3.6

2019-09-23 (Date of Last Commit)

* EmptyDrops output is now included in the ZARR output
* The GTF modification step is removed from the scatter, resulting in better performance and caching
* Memory of several tasks is increased
* The ZARR output is now compulsory and the relevant input flag has been removed
* Support for loom format has been added and a new optional flag dictates if the file is created
Documentation has been updated


# 1.3.5

2019-09-09 (Date of Last Commit)

* Increase memory for CalculateCellMetrics

# 1.3.4

2019-09-04 (Date of Last Commit)

* Increase memory

# 1.3.3

2019-08-08 (Date of Last Commit)

* Release a new patch version of Optimus with an ambitious memory allocation for CalculateCellMetrics task.
* This version and all future versions have been scientifically validated on Mouse reference version mm10 (GRCm39, Gencode M21)

# 1.3.2

2019-08-06 (Date of Last Commit)

* Update Optimus patch version (#241)

# 1.3.1

2019-07-18 (Date of Last Commit)

* This release includes fixes to the Optimus gene id outputs

# 1.3.0

2019-06-19 (Date of Last Commit)

* This release is a backwards incompatible change for the Optimus pipeline. The output matrices now contain gencode v27 gene ids in addition to gene names to comply with outputs from other pipelines and expectations of downstream services

# 1.2.0

2019-06-04 (Date of Last Commit)

* No change note

# 1.1.0

2019-05-07 (Date of Last Commit)

* No change note

# 1.0.0

2019-03-27 (Date of Last Commit)

* The first major version release for the Optimus pipeline.
