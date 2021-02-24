# 2.2.1

2020-12-07 (Date of Last Commit)

* Added library, species and organ metadata to SmartSeq2 pipeline merged loom file
* Updated the docker in LoomUtils.wdl task to 0.0.6. Updated merge_loom.py in the docker

# 2.2.0

2020-12-04 (Date of Last Commit)

* Added Gene as row attribute of the loom file
* Added CellID as column attributes of the loom file
* Updated the docker in LoomUtils.wdl task to 0.0.5. Updated create_loom_ss2.py in the docker.

# 2.1.5

2020-11-24 (Date of Last Commit)

* Made CPU, memory, and disk optional parameters for all tasks

# 2.1.4

2020-11-05 (Date of Last Commit)

* Added input checking code into HISAT tasks (called by the single sample workflow) in order to reduce number of DRS lookups

# 2.1.3

2020-10-26 (Date of Last Commit)

* Changed the SS2 single sample global attributes "input_id" and "input_name" to column attributes
* Updated the docker in LoomUtils.wdl task to 0.0.4-ss2-loom-fix-1

# 2.1.2

2020-10-13 (Date of Last Commit)

* Fixed a bug in the loom file generation script that appeared when using optional input `input_id_metadata_field`
* Updated the docker in LoomUtils.wdl task to v0.0.3

# 2.1.1

2020-10-01 (Date of Last Commit)

* Added checks for compressed fastq input files

# 2.1.0

2020-08-10 (Date of Last Commit)

### Non-breaking changes
* Added batch_name as an optional input for user provided biomaterial id
* Passed pipeline_version to output loom file  
* Added input_id_metadata_field and input_name_metadata_field as optional input

# 2.0.1

2020-07-20 (Date of Last Commit)

* Changed the imports to relative imports to support Dockstore->Terra release

# 2.0.0

2020-06-04 (Date of Last Commit)

* Removed zarr output and made loom output as default

* Loom file attribute names have changed: CellID: cell_names, Gene: gene_names and Accession: ensembl_ids

* Loom file name has changed from out.loom to "plateid".loom

* Added the expected counts in addition to TPMs in the loom matrix. 

# 1.1.0

2020-05-07 (Date of Last Commit)

* Added estimated count matrix to zarr output

# 1.0.0

2019-12-16 (Date of Last Commit)

* This is the first release of the Smart-seq2 Multi Sample workflow
