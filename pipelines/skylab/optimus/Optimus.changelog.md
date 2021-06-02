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
