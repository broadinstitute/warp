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
