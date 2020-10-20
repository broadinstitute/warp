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
