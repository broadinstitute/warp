# 1.16.0
2022-09-07 (Date of Last Commit)

* Updated task IngestOutputsToTDR in InternalTasks.wdl with new docker tag to accommodate changes for BroadInternalArrays pipeline. Change has no effect on this pipeline.

# 1.15.9
2022-06-15 (Date of Last Commit)

* Added bead_pool_manifest file to the output section of ValidateChip.wdl

# 1.15.8
2022-05-19 (Date of Last Commit)

* Patch security vulnerability in arrays-picard-private docker image
* Update arrays internal tasks, this update has no effect on this pipeline

# 1.15.7
2022-04-19 (Date of Last Commit)

* Updated to Picard version 2.26.10 and GATK version 4.2.6.1 to address log4j vulnerabilities

# 1.15.6
2022-04-14 (Date of Last Commit)

* Update base image for picard-private docker image
* Add gsutil to PATH in picard-private docker image

# 1.15.5
2022-02-025 (Date of Last Commit)

* Update to Picard 2.26.11
  * Address obscure bug in GtcToVcf -> VcfToAdpc (some variant metrics, calculated as infinite, were rendered incorrectly in the VCF)

# 1.15.4
2022-01-19  (Date of Last Commit)

* Update version of gatk in used in ValidateVariants to 4.2.4.1 (updated to log4j 2.17.1)

# 1.15.3
2022-01-18  (Date of Last Commit)

* Increase Boot disk for GATK tasks to avoid an out of disk space error

# 1.15.2
2022-01-14 (Date of Last Commit)

* Fix issue with escaping of strings/filenames with spaces embedded that occurred on older (< 57) versions of Cromwell
* Task wdls used by Validate chip were updated with changes that don't affect ValidateChip wdl

# 1.15.1
2022-01-11

* Updated picard and picard-related tasks to Picard 2.26.10
    * Address log4shell security issue (updated to log4j 2.17.1)

# 1.15.0
2021-11-17

* Updated to Picard 2.26.4
    * Changed GtcToVcf to account for zeroed-out SNPs in the calculation of GTC Call Rate. Previously the GTC Call Rate (which is stored in the VCF header) had been copied directly from the Illumina GTC File. However Illumina's calculation of the GTC Call Rate does not account for (ignore) zeroed-out SNPs, so we recalculate the GTC Call Rate, ignoring zeroed-out SNPs and use this.
    * Fixed a bug in GtcToVcf where 'SOURCE' fields read from the Illumina manifest that contain a semicolon may be incorrectly populated in the INFO field of the VCF.
* Lowered call rate threshold used in autocall to determine if a gender call can be made in order to compensate for Illumina's GTC Call Rate not accounting for zeroed-out SNPs.

# 1.14.1
2021-11-10

* Added Xmx flag (maximum heap size) to all tasks with java commands

# 1.14.0
2021-10-07

* Updated Illumina IAAP Autocall docker image to address critical vulnerability
* Changed arrays-picard-private hash and pull from correct artifactory
* Task wdls used by Validate chip were updated with changes that don't affect ValidateChip wdl
* Change outputs of ValidateChip pipeline to use python_file_naming_convention instead of CamelCase

# 1.13.3
2021-09-09

* Set the volatile=true flag for several internal tasks so they will not use call-caching

# 1.13.2
2021-08-10

* Updated GtcToVcf task in IlluminaGenotypingArrayTasks to escape fingerprint file names.
* Updated Illumina IAAP Autocall to alpine base image

# 1.13.1
2021-07-23

* Task wdls used by Validate chip were updated with changes that don't affect ValidateChip wdl

# 1.13.0
2021-05-19

* Update to use publicly released version of CreateExtendedIlluminaManifest (in Picard 2.25.5)
* Update version of Picard to 2.25.5 in order to allow GtcToVcf to support new enums in that buid (updated all picard tools to use this version)

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
