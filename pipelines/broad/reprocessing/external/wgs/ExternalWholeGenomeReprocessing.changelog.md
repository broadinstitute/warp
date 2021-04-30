# 1.3.5
2021-03-17

* Promoted VariantCalling to be a top-level workflow

# 1.3.4
2021-02-22

* Added SORTING_COLLECTION_SIZE_RATIO as an optional task input to MarkDuplicates

# 1.3.3
2021-02-08

* Calculate java memory value from the optional memory input value for CramToUnmappedBams java tasks

# 1.3.2
2021-02-02

* Minor changes to support CramToUnmappedBams as an independent versioned pipeline
    * Changed path of the relative import
    * Added 'base_file_name' as an input to CramToUnmappedBams

# 1.3.1
2020-12-21

* Passed an input bam index to several subworkflows, so the pipeline passes on singularity for sharded BQSR

# 1.3.0
2020-12-16

* Fixed error in relative import statement in Alignment subworkflow.
* Fixed syntax bug in Alignment task SamToFastqAndBwaMemAndMba

# 1.2.0
2020-10-20

* Updated GATK docker image for all tasks to [GATK 4.1.8.0](https://github.com/broadinstitute/gatk/releases/tag/4.1.8.0).
    * Numerous bug fixes and improvements
* Updated Picard docker image for all tasks to [2.23.8](https://github.com/broadinstitute/picard/releases/tag/2.23.8).
* Updated samtools to version [1.11](https://github.com/samtools/samtools/releases/tag/1.11).  Primarily for improved compression of cram files.

# 1.1.1
2020-10-01

* Removed extra trailing slash in ouput directory from cloud to cloud copy job

# 1.1.0
2020-08-18

* Added a meta section with 'allowNestedInputs' set to 'true' to allow the workflow to use with nested inputs with Cromwell 52

# 1.0.2
2020-07-31

* Update various tasks to stop using phusion/baseimage:latest docker image (it has been removed).  Start using a Google-hosted base image in it's stead.

# 1.0.1
2020-07-15

* Remove GetBWAVersion as a task and moved it to SamToFastqAndBwaMemAndMba

# 1.0
2020-04-30

* Initial release of the ExternalWholeGenomeReprocessing pipeline
