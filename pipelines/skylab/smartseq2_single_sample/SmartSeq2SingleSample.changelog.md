# 5.1.16
2023-01-23 (Date of Last Commit)

* Added "Disk" to task runtime sections to support running on Azure
* Addressed mb/gb memory specification inconsistencies in LoomUtils and CheckInput

# 5.1.15
2022-09-13 (Date of Last Commit)

* Update RSEM.wdl to use an updated RSEM docker image. This change does not affect the SmartSeq2SingleSample pipeline.

# 5.1.14
2022-09-12 (Date of Last Commit)

* Update HISAT2.wdl to use an updated HISAT2 docker image. This change does not affect the SmartSeq2SingleSample pipeline.

# 5.1.13
2022-08-23 (Date of Last Commit)

* Remove an unused script in pytools docker image.

# 5.1.12
2022-08-16 (Date of Last Commit)

* Updated LoomUtils.wdl to use a consolidated python utilities docker image. This change does not affect the SmartSeq2SingleSample pipeline.

# 5.1.11
2022-06-21 (Date of Last Commit)

* Updated main workflow name from SmartSeq2SingleCell to SmartSeq2SingleSample in the SS2 single sample pipeline. This allows the pipeline to run in the updated scala tests.

# 5.1.10
2022-04-22 (Date of Last Commit)

* Updated LoomUtils.wdl for a task in the Optimus pipeline. This change does not affect the SmartSeq2SingleSample pipeline.

# 5.1.9
2022-04-14 (Date of Last Commit)

* Updated to Picard version 2.26.10 and GATK version 4.2.6.1 to address log4j vulnerabilities
    * Two new metrics added to insert size metrics: 
        * MODE_INSERT_SIZE
        * WIDTH_OF_95_PERCENT

# 5.1.8
2022-02-25 (Date of Last Commit)

* Updated LoomUtils.wdl for a task in the Optimus pipeline. This change does not affect the SmartSeq2SingleSample pipeline.
# 5.1.7
2022-01-07 (Date of Last Commit)

* Updated LoomUtils.wdl to fix a missing metadata issue in Single Nucleus SmartSeq2 pipeline
# 5.1.6
2021-11-10 (Date of Last Commit)

* Added Xmx flag (maximum heap size) to all tasks with java commands
# 5.1.3
2022-01-07 (Date of Last Commit)

* Updated LoomUtils.wdl to fix a missing metadata issue in Single Nucleus SmartSeq2 pipeline
# 5.1.5
2021-09-13 (Date of Last Commit)

* Updated Picard.wdl and LoomUtils.wdl for Single Nucleus SmartSeq2. These changes do not affect SmartSeq2SingleSample

# 5.1.4
2021-09-02 (Date of Last Commit)

* Removed a redundant task in Picard.wdl that was use in the previous  Optimus pipeline. However, 
  that wdl also contains other Picard task that are used in the smartseq2 single sample. Therefore, 
  the smartseq2 single sample is not expected to change. 

# 5.1.3
2021-07-19 (Date of Last Commit)

* Updated SmartSeq2 to accommodate spaces in input_name

# 5.1.2

2020-05-24 (Date of Last Commit)

* Updated the LoomUtils.wdl and Picard.wdl. The changes are made for smartseq2 single nuclei and  multisample smartseq2 single nuclie pipelines, but does not impact SmartSeq2 single sample pipeline

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


