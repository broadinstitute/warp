# 2.6.18
2023-01-13 (Date of Last Commit)

* Updated remaining uses of GATK to version 4.3.0.0

# 2.6.17
2022-11-29 (Date of Last Commit)

* Updated InternalArraysTasks.wdl for BroadInternalArrays - adding `lab_batch` as workflow input.
* Updated to GATK version 4.3.0.0

# 2.6.16
2022-11-08 (Date of Last Commit)

* Removed task MakeOptionalOutputBam in Utilities.wdl, this update has no effect on this pipeline.
* Updated task IngestOutputsToTDR in InternalTasks.wdl with new optional input variable. This update has no effect on this pipeline.
* Updated task FormatArraysOutputs in InternalArrraysTasks.wdl with new docker tag to accommodate changes for BroadInternalArrays pipeline. Change has no effect on this pipeline.
* Updated task IngestOutputsToTDR in InternalTasks.wdl with new docker tag to accommodate changes for BroadInternalArrays pipeline. Change has no effect on this pipeline.

* New GCR image tag in task FormatImputationOutputs and task FormatImputationWideOutputs and task IngestOutputsToTDR.
* Removed workspace_bucket workflow parameter from BroadInternalArrays and BroadInternalImputation.

# 2.6.15
2022-09-30 (Date of Last Commit)

* Updated Picard-Python Docker image in Utilities.wdl to fix vulnerabilities.
* Updated tasks FormatArraysOutputs and IngestOutputsToTDR with GCR images instead of Dockerhub.

# 2.6.14
2022-09-07 (Date of Last Commit)

* Updated task MakeOptionalOutputBam in Utilities.wdl, this update has no effect on this pipeline
* Updated task FormatArraysOutputs in InternalArrraysTasks.wdl with new docker tag to accommodate changes for BroadInternalArrays pipeline. Change has no effect on this pipeline.
* Updated task IngestOutputsToTDR in InternalTasks.wdl with new docker tag to accommodate changes for BroadInternalArrays pipeline. Change has no effect on this pipeline.

# 2.6.13
2022-06-21 (Date of Last Commit)

* Changed QC.CheckFingerprint to QC.CheckFingerprintTask to avoid a naming conflict in the update scala tests, no effect on this pipeline

# 2.6.12
2022-06-01 (Date of Last Commit)

* Renamed the CompareVCFs task in VerifyIlluminaGenotypingArray.wdl to CompareVcfsAllowingQualityDifferences, this update has no effect on this pipeline

# 2.6.11
2022-05-19 (Date of Last Commit)

* Patch security vulnerability in arrays-picard-private docker image
* Update arrays internal tasks, this update has no effect on this pipeline

# 2.6.10
2022-05-05 (Date of Last Commit)
* Update Arrays.wdl to output bead_pool_manifest_file so it can be used in new test wdls

# 2.6.9
2022-04-19 (Date of Last Commit)
* Updated to GATK version 4.2.6.1 to address log4j vulnerabilities

# 2.6.8
2022-04-04 (Date of Last Commit)

* Update task CopyFilesFromCloudToCloud in Utilities.wdl, this update has no effect on this pipeline

# 2.6.7
2022-03-24 (Date of Last Commit)

* Update base image for picard-private docker image
* Add gsutil to PATH in picard-private docker image

# 2.6.6
2022-02-25 (Date of Last Commit)

* Update to Picard 2.26.11
    * Address obscure bug in GtcToVcf -> VcfToAdpc (some variant metrics, calculated as infinite, were rendered incorrectly in the VCF)

# 2.6.5
2022-02-01  (Date of Last Commit)

* Address memory usage in CheckFingerprint task to allow anough headroom for the VM

# 2.6.4
2022-01-19  (Date of Last Commit)

* Update version of gatk in used in tasks in IlluminaGenotypingArray WDL to 4.2.4.1 (updated to log4j 2.17.1)

# 2.6.3
2022-01-18  (Date of Last Commit)

* Increase Boot disk for GATK tasks to avoid an out of disk space error

# 2.6.2
2022-01-14  (Date of Last Commit)

* Fix issue with escaping of strings/filenames with spaces embedded that occurred on older (< 57) versions of Cromwell
* Refactor to move CheckFingerprint functionality into separate task

# 2.6.1
2022-01-11

* Updated picard and picard-related tasks to Picard 2.26.10
    * Address log4shell security issue (updated to log4j 2.17.1)

# 2.6.0
2021-11-17

* Updated to Picard 2.26.4
    * Changed GtcToVcf to account for zeroed-out SNPs in the calculation of GTC Call Rate. Previously the GTC Call Rate (which is stored in the VCF header) had been copied directly from the Illumina GTC File. However Illumina's calculation of the GTC Call Rate does not account for (ignore) zeroed-out SNPs, so we recalculate the GTC Call Rate, ignoring zeroed-out SNPs and use this.
    * Fixed a bug in GtcToVcf where 'SOURCE' fields read from the Illumina manifest that contain a semicolon may be incorrectly populated in the INFO field of the VCF.
* Lowered call rate threshold used in autocall to determine if a gender call can be made in order to compensate for Illumina's GTC Call Rate not accounting for zeroed-out SNPs.

# 2.5.3
2021-11-15

* Task wdls used by the Arrays pipeline were updated with changes that don't affect Arrays wdl

# 2.5.2
2021-11-10

* Added Xmx flag (maximum heap size) to all tasks with java commands

# 2.5.1
2021-10-25

* Make fingerprint retrieval and storage tasks use max_retries to enable recovery from transient failures
* Modified Arrays pipeline to not read fingerprints for control samples from the Mercury Fingerprint Store.

# 2.5.0
2021-10-07

* Enabled pipeline to lookup the extended_illumina_manifest_file using an alternate method
    * If the path to the file is not provided, it will look in the arrays_metadata_path for a map file that contains a mapping of chip to extended_illumina_manifest
* Enabled pipeline to lookup the cluster_file using an alternate method (using the arrays_metadata_path and cluster_filename)
* Enabled pipeline to lookup the gender_cluster_file using an alternate method (using the arrays_metdata_path and gender_cluster_filename)
* Enabled pipeline to lookup the zcall_thresholds_file using an alternate method (using the arrays_metdata_path and zcall_thresholds_filename)
* Enabled pipeline to lookup the genotype control data using an alternate method (using the arrays_control_data_path and control_sample_name)
* Modified pipeline to NOT write fingerprints for control samples to the Mercury Fingerprint Store.
* Change outputs of Arrays and pipeline to use python_file_naming_convention instead of CamelCase
* Removed the volatile=true flag from UploadFingerprintToMercury

# 2.4.2
2021-09-22

* Enabled pipeline to lookup the bead_pool_manifest_file using an alternate method (using the arrays_metadata_path and bead_pool_manifest_filename)

# 2.4.1
2021-09-09

* Changed default threshold for passing control (HapMap) genotype concordance from 0.98 to 0.95
* Modified pipeline to automatically generate the analysis_version_number if it is not supplied as an input.
* Modified pipeline to make several inputs optional:
    * sample_id
    * participant_id
    * collaborator_participant_id
    * lab_batch
    * product_family
    * product_name
    * product_order_id
    * product_part_number
* Set the volatile=true flag for several internal tasks so they will not use call-caching

# 2.4.1
2021-08-25
* Updated arrays-picard-private docker image to address critical vulnerability
* Changed arrays-picard-private hash and pull from correct artifactory

# 2.4.0
2021-08-05

* Enable pipeline to (optionally) pull and push fingerprints from/to the Mercury Fingerprint Store

# 2.3.6
2021-08-02

* Increased the version number to make new release tag for Dockstore 

# 2.3.5
2021-07-29

* Updated documentation to describe changes to inputs and outputs

# 2.3.4
2021-07-28

* Set a default value for product_type so that it can be safely omitted from input file

# 2.3.3
2021-07-22

* Have pipeline take the values supplied in 'params.txt' input file as optional top-level inputs. First step towards removal
* Provide params.txt file as output of pipeline.
* Set default call rate threshold of pipeline to 0.98

# 2.3.2
2021-7-19

* Updated Illumina IAAP Autocall to alpine base image
* Make chip_well_barcode and analysis_version_number available as outputs of the WDL.
  
# 2.3.1
2021-05-19

* Update version of Picard to 2.25.5 in order to allow GtcToVcf (used in IlluminaGenotypingArray subworkflow) to support new enums in that buid (updated all picard tools to use this version)

# 2.3.0
2020-10-07

* Added use of BafRegress to the pipeline.  BafRegress detects and estimates sample contamination using B allele frequency data from Illumina genotyping arrays using a regression model.
* Updated all internal tasks to use the latest version of picard-private as best practice.

# 2.2.0
2020-10-01

* Updated task definitions to include a new tool not currently used in Arrays wdl

# 2.1.0
2020-08-18

* Added a meta section with 'allowNestedInputs' set to 'true' to allow the workflow to use with nested inputs with Cromwell 52

# 2.0.1
2020-08-31

* Fixed regression in the fix for VcfToMercuryFingerprintJson task

# 2.0
2020-07-31

* Fix bug in VcfToMercuryFingerprintJson task where the metrics file is being parsed incorrectly
  if there is a space in the sample name value 
  
# 1.9.1
2020-07-31

* Have md5sum task stop using phusion/baseimage:latest docker image (it has been removed).  Start using a Google-hosted base image in it's stead.

# 1.9
2020-06-23

* Updated to use new version of GtcToVcf that uses bpm, rather than bpm.csv file.
* Update database-accessing tasks to no longer use cloud proxy.

# 1.8
2020-06-10

* Updated Picard-Private to support population of additional column 'PRODUCT_TYPE' in CHIP_WELL_BARCODE_INDEX table by the UpdateChipWellBarcodeIndex CLP
* Updated Picard-Private to support population of additional column 'PIPELINE_VERSION' in ARRAYS_QC table by the UploadArraysMetrics CLP

# 1.7

2020-05-15

*Added md5s for red idats, green idats and vcf input/outputs of the arrays pipeline.

# 1.6
2020-03-05

* Allow use of optional 'contamination_controls_vcf' which lets the contamination checking tool 'VerifyIDIntensity' be run in multi-sample mode, improving contamination estimation.  
* Updated Picard to 2.22.0
    * Support multiple sample VCFs in VcfToAdpc tool

# 1.5

2019-12-19

* Updated Picard to 2.21.6
    * Updated CreateVerifyIDIntensityContaminationMetricsFile to handle negative LLK values
    * Corrected float output of GtcToVcf to avoid downstream problems in vcf parsing 
* Improve accuracy of Genotype Concordance with control samples by excluding filtered sites from VCF prior to calling GenotypeConcordance tool.

# 1.4

2019-12-16

* Adjusted memory parameters to avoid Google's new e2 instances because there are not enough machines to satisfy our production use case.

# 1.3

2019-11-21

* Renamed file 'ArraysWf.wdl' to 'Arrays.wdl'

# 1.2
* Moved two docker images used in the IlluminaGenotypingArray pipeline from private registry to public

# 1.1

## Bug fix, updated documentation
* Update WDL tasks to use latest Picard docker image
    * Bug fix in VcfToAdpc (Allow for null values of Normalized X and Y intensity)  https://github.com/broadinstitute/picard/pull/1415
* Bug fix in determining autocall gender.  For calling IAAP/gencall with a gender cluster file, we now use an autosomal call rate threshold < 0. This is because IAAP does not calculate gender if the autosomal call rate is below some threshold.  But for Broad gender cluster files, the lab has zeroed out the autosomal SNPs and the autosomal call rate will be 0.
* Update documentation

# 1.0
Initial public release of the Single-sample Arrays pipeline

## Detailed Release Notes

* We have changed from a version of Illumina Autocall shared under NDA to a public Illumina tool called [IAAP](https://support.illumina.com/downloads/iaap-genotyping-cli.html). This resulted in no changes to genotype calls in our tests.
    * The reported version number in the GTC will change.
    * IAAP requires more evidence for gender calling and will report “Unknown” more often than Autocall did.
    * Metrics in the GTC file from IAAP are different from Autocall metrics.
* For chips that are called with zCall, we have removed a filter applied to sites where the zCall-called genotype differs from that of the autocall genotype. Differences can still be found by using other fields in the VCF. This will result in approximately 5% more unfiltered calls on chips that are zCalled (chip dependent).
* Contamination - [VerifyIDIntensity](https://github.com/gjun/verifyIDintensity) is an externally developed tool that we have added to the pipeline. The pipeline outputs a new metrics file regarding the contamination of the sample using the outputs of this tool. 
* Metrics - There are 2 new metrics and renaming of 4 metrics.
    * New
        * `NUM_ASSAYS` - total num of assays on the VCF
        * `NUM_ZEROED_OUT_ASSAYS` - number of zeroed out assays (these are assays that are marked as unusable in the cluster file)
    * Renamed
        * `FILTERED_SNPS' -> 'NUM_FILTERED_ASSAYS`
        * `TOTAL_ASSAYS' -> 'NUM_NON_FILTERED_ASSAYS`
        * `TOTAL_SNPS' -> 'NUM_SNPS`
        * `TOTAL_INDELS' -> 'NUM_INDELS`

The Data Sciences Platform is dedicated to Open Source software and open science. To support this pipeline and that mission, we have Open Sourced the following tools in Picard: 

* [GtcToVcf](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_arrays_GtcToVcf.php) - Class to convert a GTC file and a BPM file to a VCF file.
* [MergePedIntoVcf](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_arrays_MergePedIntoVcf.php) - Class to take genotype calls from a ped file output from zCall and merge them into a VCF generated by autocall.
* [CollectArraysVariantCallingMetrics](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_arrays_CollectArraysVariantCallingMetrics.php) - Collects summary and per-sample metrics about variant calls in a VCF file.
* [VcfToAdpc](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_arrays_VcfToAdpc.php) - A simple program to convert a Genotyping Arrays VCF to an ADPC file (Illumina intensity data file).
* [CreateVerifyIdIntensityContaminationMetricsFile](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_arrays_CreateVerifyIDIntensityContaminationMetricsFile.php) - A simple program to create a standard picard metrics file from the output of VerifyIDIntensity
