# 1.12.0
2020-10-07

* Updated task definitions to include a new tool not currently used in ValidateChip wdl
* Updated all internal tasks to use the latest version of picard-private as best practice.

# 1.11.0
2020-10-01

* Updated task definitions to include a new tool not currently used in ValidateChip wdl

# 1.10.0
2020-08-18

* Added a meta section with 'allowNestedInputs' set to 'true' to allow the workflow to use with nested inputs with Cromwell 52

# 1.9.1
2020-08-31

* Update name of extended Illumina manifest file to be commensurate with version used in picard-private.
# 1.9
2020-07-31

* Fix bug in GenotypeConcordance task where the metrics file is being parsed incorrectly
  if there is a space in the sample name value 

# 1.8.1
2020-07-31

* Update various tasks to stop using phusion/baseimage:latest docker image (it has been removed).  Start using a Google-hosted base image in it's stead.

# 1.8
2020-06-23

* Updated to use new version of GtcToVcf that uses bpm, rather than bpm.csv file
* Updated Picard-Private for internal arrays tasks.

# 1.7
2020-06-10

* Updated Picard-Private for use in task not used in ValidateChip

# 1.6
2020-05-15

*Changes to arrays pipeline to add md5sums of the idats and vcf

# 1.5
2020-03-05

* Updated Picard to 2.22.0
    * Support multiple sample VCFs in VcfToAdpc tool

# 1.4

2019-12-19

* Updated Picard to 2.21.6
    * Updated CreateVerifyIDIntensityContaminationMetricsFile to handle negative LLK values
    * Corrected float output of GtcToVcf to avoid downstream problems in vcf parsing 
* Refactor to use common SelectVariants task in this and other pipeline wdls

# 1.3

2019-12-16

* Adjusted memory parameters to avoid Google's new e2 instances because there are not enough machines to satisfy our production use case.

# 1.2
Moved two docker images used by ValidateChip from the private registry to public

# 1.1

## Updated picard dockers, updated documentation
* Update WDL tasks to use latest Picard docker image
* Update documentation

# 1.0
Initial release of the ValidateChip pipeline
