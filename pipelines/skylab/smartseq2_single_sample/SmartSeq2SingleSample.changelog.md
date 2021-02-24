# 5.1.1

2020-12-07 (Date of Last Commit)

* Updated the docker in LoomUtils.wdl task to 0.0.6. Updated merge_loom.py in the docker

# 5.1.0

2020-12-04 (Date of Last Commit)

* Added Gene as row attribute of the loom file
* Added CellID as column attributes of the loom file
* Updated the docker in LoomUtils.wdl task to 0.0.5. Updated create_loom_ss2.py script in the docker.

# 5.0.5

2020-11-24 (Date of Last Commit)

* Made CPU, memory, and disk optional parameters for all tasks

# 5.0.4

2020-11-05 (Date of Last Commit)

* Pushed input checking code down into HISAT tasks in order to reduce number of DRS lookups

# 5.0.3

2020-10-26 (Date of Last Commit)

* Changed the SS2 single sample global attributes "input_id" and "input_name" to column attributes
* Updated the docker in LoomUtils.wdl task to 0.0.4-ss2-loom-fix-1

# 5.0.2

2020-10-13 (Date of Last Commit)

* Fixed a bug in the loom file generation script that appeared when using optional input `input_id_metadata_field`
* Updated the docker in LoomUtils.wdl task to v0.0.3

# 5.0.1

2020-10-01 (Date of Last Commit)

* Added check to see if input fastq files are compressed in HISAT2.wdl task

# 5.0.0

2020-08-10 (Date of Last Commit)
### Breaking changes 
* Changed sample_name to input_id
### Non-breaking changes
* Added input_name as an optional input for user provided biomaterial id.
* Passed pipeline_version to output loom file
* Added input_id_metadata_field and input_name_metadata_field as optional input


# 4.0.1

2020-07-20 (Date of Last Commit)

* Changed the imports to relative imports to support Dockstore->Terra release

# 4.0.0

2020-06-04 (Date of Last Commit)

* Added loom output and removed the zarr output

* Loom file attribute names have changed: CellID: cell_names, Gene: gene_names and Accession: ensembl_ids

* Loom file name has changed from out.loom to "sample_id".loom

* Added expected counts in addition to the TPMs in the loom matrix

# 3.1.0

2020-05-07 (Date of Last Commit)

* Added estimated count matrix to zarr output

# 3.0.0

2019-12-16 (Date of Last Commit)

* Added single end support for the Smart-seq2 Single Sample pipeline
* Removed a previous version of a standalone single-end pipeline


# 2.5.0

2019-11-07 (Date of Last Commit)

* Removed max_retries parameter so that a default value can be set by the workflow options

# 2.4.0

2019-05-24 (Date of Last Commit)

* No release note available

# 2.3.0

2019-05-07 (Date of Last Commit)

* No release note available

# 2.2.0

2019-01-22 (Date of Last Commit)

* No release note available

# 2.1.0

2018-10-16 (Date of Last Commit)

* No release note available

# 2.0.0 

2018-10-09 (Date of Last Commit)

* No release note available

# 1.0.0 

2018-06-14 (Date of Last Commit)

* No release note available


