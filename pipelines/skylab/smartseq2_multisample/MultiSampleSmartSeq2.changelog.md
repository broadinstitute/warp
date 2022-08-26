# 2.2.14
2022-08-23 (Date of Last Commit)

* Remove an unused script in pytools docker image.

# 2.2.13
2022-08-16 (Date of Last Commit)

* Update LoomUtils.wdl to use updated docker images. This change does not affect the MultiSampleSmartSeq2 pipeline.

# 2.2.12
2022-06-22 (Date of Last Commit)

* Updated main workflow name from SmartSeq2SingleCell to SmartSeq2SingleSample in the SS2 single sample pipeline. This allows the pipeline to run in the updated scala tests.

# 2.2.11
2022-04-22 (Date of Last Commit)

* Updated LoomUtils.wdl for a task in the Optimus pipeline. This change does not affect the MultiSampleSmartSeq2 pipeline.

# 2.2.10
2022-04-14 (Date of Last Commit)

* Updated to Picard version 2.26.10 and GATK version 4.2.6.1 to address log4j vulnerabilities

# 2.2.9
2022-02-25 (Date of Last Commit)

* Updated LoomUtils.wdl for a task in the Optimus pipeline. This change does not affect the MultiSampleSmartSeq2 pipeline. 

# 2.2.8
2022-01-07 (Date of Last Commit)

* Updated LoomUtils.wdl to fix a missing metadata issue in Single Nucleus SmartSeq2 pipeline

# 2.2.7
2021-11-10 (Date of Last Commit)

* Added Xmx flag (maximum heap size) to all tasks with java commands

# 2.2.6
2021-09-13 (Date of Last Commit)

* Updated Picard.wdl and LoomUtils.wdl for Single Nucleus SmartSeq2. These changes do not affect MultiSampleSmartSeq2

# 2.2.5
2021-09-02 (Date of Last Commit)

* Removed a redundant task in Picard.wdl that was use in the previous  Optimus pipeline. However, 
  that wdl also contains other Picard task that are used in the smartseq2 single sample. Therefore, 
  the smartseq2 single sample is not expected to change. 

# 2.2.4
2021-08-02 (Date of Last Commit)

* Increased the version number to make new release tag for Dockstore 

# 2.2.3
2021-07-19 (Date of Last Commit)

* Updated SmartSeq2 to accommodate spaces in input_name

# 2.2.2

2020-05-24 (Date of Last Commit)

* Added a task to Picard.wdl
* Updated the docker in LoomUtils.wdl task to 0.0.7. 

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
