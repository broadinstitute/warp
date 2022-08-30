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
