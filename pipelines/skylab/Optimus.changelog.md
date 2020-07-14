# optimus_v3.0.0
2020-06-10 (Date of Last Commit)

* Removed the Zarr formatted matrix and metrics outputs and replaced with Loom
* Removed EmptyDrops for sn_rna mode
* Updated the Loom file attribute names: CellID to cell_names, Gene to gene_names, and Accession to ensembl_ids
* Added metrics for mitochondrial reads
* Added an optional input for the BAM basename; this input is listed as ‘bam_output_basename’and the default is 'sample_id'

# optimus_v2.0.0
2020-02-08 (Date of Last Commit)

* Fixed a bug that resulted in emptyDrops output being incorrect
* Updated the workflow to WDL 1.0

# optimus_v1.4.0

2019-11-08 (Date of Last Commit)

* Addition of support for V3 chemistry
* Addition of input parameter validation step
* Greatly improved documentation
* Improvements to ZARR output

# optimus_v1.3.6

2019-09-23 (Date of Last Commit)

* EmptyDrops output is now included in the ZARR output
* The GTF modification step is removed from the scatter, resulting in better performance and caching
* Memory of several tasks is increased
* The ZARR output is now compulsory and the relevant input flag has been removed
* Support for loom format has been added and a new optional flag dictates if the file is created
Documentation has been updated


# optimus_v1.3.5

2019-09-09 (Date of Last Commit)

* Increase memory for CalculateCellMetrics

# optimus_v1.3.4

2019-09-04 (Date of Last Commit)

* Increase memory

# optimus_v1.3.3

2019-08-08 (Date of Last Commit)

* Release a new patch version of Optimus with an ambitious memory allocation for CalculateCellMetrics task.
* This version and all future versions have been scientifically validated on Mouse reference version mm10 (GRCm39, Gencode M21)

# optimus_v1.3.2

2019-08-06 (Date of Last Commit)

* Update Optimus patch version (#241)

# optimus_v1.3.1

2019-07-18 (Date of Last Commit)

* This release includes fixes to the Optimus gene id outputs

# optimus_v1.3.0

2019-06-19 (Date of Last Commit)

* This release is a backwards incompatible change for the Optimus pipeline. The output matrices now contain gencode v27 gene ids in addition to gene names to comply with outputs from other pipelines and expectations of downstream services

# optimus_v1.2.0

2019-06-04 (Date of Last Commit)

* No change note

# optimus_v1.1.0

2019-05-07 (Date of Last Commit)

* No change note

# optimus_v1.0.0

2019-03-27 (Date of Last Commit)

* The first major version release for the Optimus pipeline.


