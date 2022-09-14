# 5.5.4
2022-09-01 (Date of Last Commit)

* Update CheckInputs.wdl to use a lightweight alpine-bash image.

# 5.5.3
2022-08-23 (Date of Last Commit)

* Remove an unused script in pytools docker image and removed unused ConvertStarOutputs task.

# 5.5.2
2022-08-16 (Date of Last Commit)

* Updated LoomUtils.wdl and StarAlign.wdl to use an updated python utilities docker.

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
