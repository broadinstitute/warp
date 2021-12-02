# 1.0.4
2021-11-22 (Date of Last Commit)

* Removed haplotype_database as a pipeline input since it is not used

# 1.0.3
2021-11-19 (Date of Last Commit)

* Security patch to bcftools-vcftools docker image used in Imputation pipeline
# 1.0.2
2021-11-15 (Date of Last Commit)

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
