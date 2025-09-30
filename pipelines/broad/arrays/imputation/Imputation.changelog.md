# 1.1.22
2025-09-18 (Date of Last Commit)

* Update CountSamples task to use SSD for more efficient localization of large files on Google Batch.

# 1.1.21
2025-09-03 (Date of Last Commit)

* Update UpdateHeader task to add an optional pipeline_header_line input that when supplied, will add a header line containing 
this value to the header of the output vcf.  Currently not used by this pipeline

# 1.1.20
2025-08-26 (Date of Last Commit)

* Update tasks to use maxRetries of 1 instead of 2. 1 retry is sufficient for transient errors and helps reduce costs.

# 1.1.19
2025-07-31 (Date of Last Commit)

* Add maxRetries to tasks to help with performance on Google Batch.
* Use SSD for tasks that localize large files to help with performance on Google Batch.

# 1.1.18
2025-04-07 (Date of Last Commit)

* Update Imputation Tasks to not request external IP addresses (use runtime attribute 'noAddress: true').

# 1.1.17
2025-04-01 (Date of Last Commit)

* Update Imputation Tasks to use dockers (ubuntu and tidyverse) that have been moved to GAR.

# 1.1.16
2025-02-24 (Date of Last Commit)

* Updated runtime parameters in some ImputationTasks, and added an explicit definition of a vcf_index.

# 1.1.15
2024-11-04 (Date of Last Commit)

* Updated to GATK version 4.6.1.0

# 1.1.14
2024-09-06 (Date of Last Commit)

* Updated to GATK version 4.6.0.0

# 1.1.13
2024-05-21 (Date of Last Commit)

* Updated GermlineVariantDiscovery, BamProcessing, DragenTasks, Qc, and Utilities tasks to allow multi-cloud dockers. This change does not affect this pipeline.

# 1.1.12
2023-12-18 (Date of Last Commit)

* Updated to GATK version 4.5.0.0

# 1.1.11
2023-08-01 (Date of last Commit)

* Moved ReplaceHeader to its own scatter to remove dependency between the two nested scatters to help with wall clock time
* Updated Eagle docker to address security vulnerabilities, this has no effect on this pipeline 

# 1.1.10
2023-03-03 (Date of Last Commit)

* Adjusted disk size calculation in SplitMultiSampleVcf

# 1.1.9
2023-01-13 (Date of Last Commit)

* Updated remaining GATK uses to version 4.3.0.0

# 1.1.8
2022-12-16 (Date of Last Commit)

* Updated to GATK version 4.3.0.0

# 1.1.7
2022-12-01 (Date of Last Commit)

* Updated BCFTools/VCFTools Docker image

# 1.1.6
2022-11-10 (Date of Last Commit)

* Added meta section to allowNestedInputs, this will allow task level inputs to be set in the inputs json
* Removed task MakeOptionalOutputBam in Utilities.wdl, this update has no effect on this pipeline

# 1.1.5
2022-09-30 (Date of Last Commit)

* Updated BCFTools/VCFTools and Minimac4 Docker images to fix vulnerabilities.
* Updated tasks FormatImputationOutputs, FormatImputationWideOutputs, and IngestOutputsToTDR with GCR images instead of Dockerhub.

# 1.1.4
2022-08-23 (Date of Last Commit)

* Updated BCFTools/VCFTools docker image

# 1.1.3
2022-08-03 (Date of Last Commit)

* Updated BCFTools/VCFTools Minimac4 Docker images

# 1.1.2
2022-07-15 (Date of Last Commit)

* Updated task MakeOptionalOutputBam in Utilities.wdl, this update has no effect on this pipeline

# 1.1.1
2022-06-01 (Date of Last Commit)

* Renamed the CompareVCFs task in VerifyImputation.wdl to CompareVcfsAllowingQualityDifferences, this update has no effect on this pipeline

# 1.1.0 
2022-04-21 (Date of Last Commit)

* Update QC metrics calculation for Imputation pipeline to only evaluate sites with MAF
* Updated to Picard version 2.26.10 and GATK version 4.2.6.1 to address log4j vulnerabilities

# 1.0.9
2022-04-14 (Date of Last Commit)

* Security patch to bcftools-vcftools and minimac4 docker images

# 1.0.8
2022-04-08 (Date of Last Commit)

* Fix Imputation scaling issue by moving problematic tasks inside of a scatter call
* Update task CopyFilesFromCloudToCloud in Utilities.wdl, this update has no effect on this pipeline

# 1.0.7
2022-03-08 (Date of Last Commit)

* Add pipefail to imputation tasks to ensure that they don't fail silently

# 1.0.6
2022-03-01 (Date of Last Commit)

* Security patch to bcftools-vcftools and minimac4 docker images

# 1.0.5
2022-01-19 (Date of Last Commit)

* Security patch to bcftools-vcftools and minimac4 docker images

# 1.0.4
2021-11-22 (Date of Last Commit)

* Removed haplotype_database as a pipeline input since it is not used

# 1.0.3
2021-11-19 (Date of Last Commit)

* Security patch to bcftools-vcftools docker image used in Imputation pipeline
# 1.0.2
∑2021-11-15 (Date of Last Commit)

* Task wdls used by the Imputation pipeline were updated with changes that don't affect Imputation wdl

# 1.0.1
2021-11-10 (Date of Last Commit)
* Added Xmx flag (maximum heap size) to all tasks with java commands
* Added a step to MergeSingleSampleVcfs task to copy index files next to VCFs before running merge

    * This is needed for runnning in Broad production on data stored in TDR


# 1.0.0
2021-10-14 (Date of Last Commit)

* Initial public release of the Imputation pipeline. Read more in the [Imputation pipeline overview](https://broadinstitute.github.io/warp/docs/Pipelines/Imputation_Pipeline/README).

  * The Imputation pipeline imputes missing genotypes from either a multi-sample VCF or an array of single sample VCFs using a large genomic reference panel. It is based on the Michigan Imputation Server pipeline. Overall, the pipeline filters, phases, and performs imputation on a multi-sample VCF. It outputs the imputed VCF along with key imputation metrics.
