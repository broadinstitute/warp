# 1.12.27
2026-01-21 (Date of Last Commit)

* Moved inputs into new Google buckets. This change does not affect the outputs of the pipeline

# 1.12.26
2025-10-09 (Date of Last Commit)

* Modified the ReblockGVCF.wdl to use bash's basename instead of WDL's basename; this does not affect the outputs of this pipeline

# 1.12.25
2025-08-11 (Date of Last Commit)

* Reorganized pipelines into the wdl directory

# 1.12.24
2024-11-04 (Date of Last Commit)

* Updated to GATK version 4.6.1.0

# 1.12.23
2024-10-28 (Date of Last Commit)

* Updated the docker in the ValidateVCF task; this does not affect this pipeline

# 1.12.22
2024-09-06 (Date of Last Commit)

* Updated to GATK version 4.6.0.0

# 1.12.21
2024-08-02 (Date of Last Commit)

* The ubuntu_16_0_4 docker image version was pinned instead of using the latest tag; this does not affect the outputs of the pipeline  

# 1.12.20
2024-06-12 (Date of Last Commit)

* ValidateVcfs is more robust to larger inputs; this does not affect this pipeline

# 1.12.19
2024-07-09 (Date of Last Commit)

* Updated tasks GermlineVariantDiscovery.wdl and QC.wdl to allow multi-cloud dockers

# 1.12.18
2024-07-01 (Date of Last Commit)

* CalculateReadGroupChecksum requires more memory and disk; this does not affect this pipeline

# 1.12.17
2024-03-26 (Date of Last Commit)

* ValidateVcfs requires less memory when run without interval list; this does not affect this pipeline

# 1.12.16
2023-12-18 (Date of Last Commit)

* Updated to GATK version 4.5.0.0.

# 1.12.15
2023-12-08 (Date of Last Commit)

* ValidateVcfs now has optional memory parameter; this does not affect this pipeline

# 1.12.14
2023-11-29 (Date of Last Commit)

* ValidateVcfs can now take in VCF as calling_interval_list that is in a separate location from its index; this does not affect this pipeline

# 1.12.13
2023-03-30 (Date of Last Commit)

* CheckFingerprint can allow LOD 0

# 1.12.12
2023-01-13 (Date of Last Commit)

* Updated remaining uses of GATK to version 4.3.0.0

# 1.12.11
2022-11-09 (Date of Last Commit)

* Updated to GATK version 4.3.0.0

# 1.12.10
2022-06-21 (Date of Last Commit)

* Changed QC.CheckFingerprint to QC.CheckFingerprintTask to avoid a naming conflict in the update scala tests, no effect on this pipeline

# 1.12.9
20222-06-01  (Date of Last Commit)

* Renamed the CompareVCFs task in VerifyIlluminaGenotypingArray.wdl to CompareVcfsAllowingQualityDifferences, this update has no effect on this pipeline

# 1.12.8
20222-04-19  (Date of Last Commit)

* Updated to GATK version 4.2.6.1 to address log4j vulnerabilities

# 1.12.7
20222-04-14  (Date of Last Commit)

* Task wdls used by the IlluminaGenotypingArray pipeline were updated with changes that don't affect the IlluminaGenotypingArray pipeline itself

# 1.12.6
2022-02-025 (Date of Last Commit)

* Update to Picard 2.26.11
  * Address obscure bug in GtcToVcf -> VcfToAdpc (some variant metrics, calculated as infinite, were rendered incorrectly in the VCF)

# 1.12.5
2022-02-01 (Date of Last Commit)

* Addressed memory usage in CheckFingerprint task to allow sufficient headroom for the VM

# 1.12.4
2022-01-19  (Date of Last Commit)

* Update version of gatk in used in SubsetArrayVCF to 4.2.4.1 (updated to log4j 2.17.1)

# 1.12.3
2022-01-18  (Date of Last Commit)

* Increase Boot disk for GATK tasks to avoid an out of disk space error

# 1.12.2
2022-01-14 (Date of Last Commit)

* Fix issue with escaping of strings/filenames with spaces embedded that occurred on older (< 57) versions of Cromwell
* Refactor to move CheckFingerprint functionality into new task

# 1.12.1
2022-01-11

* Updated picard and picard-related tasks to Picard 2.26.10
    * Address log4shell security issue (updated to log4j 2.17.1)

# 1.12.0
2021-11-17

* Updated to Picard 2.26.4
    * Changed GtcToVcf to account for zeroed-out SNPs in the calculation of GTC Call Rate. Previously the GTC Call Rate (which is stored in the VCF header) had been copied directly from the Illumina GTC File. However Illumina's calculation of the GTC Call Rate does not account for (ignore) zeroed-out SNPs, so we recalculate the GTC Call Rate, ignoring zeroed-out SNPs and use this.
    * Fixed a bug in GtcToVcf where 'SOURCE' fields read from the Illumina manifest that contain a semicolon may be incorrectly populated in the INFO field of the VCF.
* Lowered call rate threshold used in autocall to determine if a gender call can be made in order to compensate for Illumina's GTC Call Rate not accounting for zeroed-out SNPs.

# 1.11.7
2021-11-10

* Added Xmx flag (maximum heap size) to all tasks with java commands

# 1.11.6
2021-10-01

* Changed the way the version of autocall is returned for the case of arrays that fail gencall
* Updated Illumina IAAP Autocall and Zcall docker images to address critical vulnerability

# 1.11.5
2021-09-08

* Changed default threshold for passing control (HapMap) genotype concordance from 0.98 to 0.95

# 1.11.4
2021-08-10

* Updated GtcToVcf task in IlluminaGenotypingArrayTasks to escape fingerprint file names.
* Updated Illumina IAAP Autocall to alpine base image

# 1.11.3
2021-08-02

* Increased the version number to make new release tag for Dockstore 

# 1.11.2
2021-07-19

* Make chip_well_barcode and analysis_version_number available as outputs of the WDL.

# 1.11.1
2021-05-19

* Update version of Picard to 2.25.5 in order to allow GtcToVcf to support new enums in that buid (updated all picard tools to use this version)

# 1.11.0
2020-10-01

* Added use of BafRegress to the pipeline.  BafRegress detects and estimates sample contamination using B allele frequency data from Illumina genotyping arrays using a regression model

# 1.10.0
2020-08-18

* Added a meta section with 'allowNestedInputs' set to 'true' to allow the workflow to use with nested inputs with Cromwell 52

# 1.9
2020-07-31

* Fixed a bug in CollectArraysVariantCallingMetrics and GenotypeConcordance where the metrics file is being parsed incorrectly if there is a space in the sample name value 
  
# 1.8.1
2020-07-31

* Removed the phusion/baseimage:latest docker image from the md5sum task and updated to a Google-hosted base image

# 1.8
2020-06-23

* Updated to use new version of GtcToVcf that uses bpm, rather than bpm.csv file

# 1.7
2020-06-10

*  Updated Picard-Private to support population of additional column 'PIPELINE_VERSION' in ARRAYS_QC table

# 1.6
2020-05-15

* Added md5sums for red idat, green idat, and vcf input/outputs of the arrays pipeline
 
# 1.5
2020-03-05

* Allow use of optional 'contamination_controls_vcf' which lets the contamination checking tool 'VerifyIDIntensity' be run in multi-sample mode, improving contamination estimation.  
* Updated Picard to 2.22.0
    * Support multiple sample VCFs in VcfToAdpc tool

# 1.4
2019-12-19

* Updated Picard to 2.21.6
    * Updated CreateVerifyIDIntensityContaminationMetricsFile to handle negative LLK values
    * Corrected float output of GtcToVcf to avoid downstream problems in vcf parsing 
* Improve accuracy of Genotype Concordance with control samples by excluding filtered sites from VCF prior to calling GenotypeConcordance tool.

# 1.3

2019-12-16

* Adjusted memory parameters to avoid Google's new e2 instances because there are not enough machines to satisfy our production use case.

# 1.2
All docker images used by the IlluminaGenotypingArray pipeline are now public

# 1.1

## Bug fix, updated documentation
* Update WDL tasks to use latest Picard docker image
    * Bug fix in VcfToAdpc (Allow for null values of Normalized X and Y intensity)  https://github.com/broadinstitute/picard/pull/1415
* Update documentation

# 1.0
Initial public release of the IlluminaGenotypingArray pipeline

## Detailed Release Notes:

The Broad Data Sciences Platform is dedicated to Open Source software and open science. To support this pipeline and that mission, we have Open Sourced the following tools in Picard: 

* [GtcToVcf](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_arrays_GtcToVcf.php) - Class to convert a GTC file and a BPM file to a VCF file.
* [MergePedIntoVcf](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_arrays_MergePedIntoVcf.php) - Class to take genotype calls from a ped file output from zCall and merge them into a VCF generated by autocall.
* [CollectArraysVariantCallingMetrics](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_arrays_CollectArraysVariantCallingMetrics.php) - Collects summary and per-sample metrics about variant calls in a VCF file.
* [VcfToAdpc](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_arrays_VcfToAdpc.php) - A simple program to convert a Genotyping Arrays VCF to an ADPC file (Illumina intensity data file).
* [CreateVerifyIdIntensityContaminationMetricsFile](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_arrays_CreateVerifyIDIntensityContaminationMetricsFile.php) - A simple program to create a standard picard metrics file from the output of VerifyIDIntensity

