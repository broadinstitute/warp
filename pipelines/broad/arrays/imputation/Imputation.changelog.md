# 1.0.1
2021-11-10

* Added Xmx flag (maximum heap size) to all tasks with java commands

# 1.0.0
2021-10-14 (Date of Last Commit)

* Initial public release of the Imputation pipeline. Read more in the [Imputation pipeline overview](https://broadinstitute.github.io/warp/docs/Pipelines/Imputation_Pipeline/README).

  * The Imputation pipeline imputes missing genotypes from either a multi-sample VCF or an array of single sample VCFs using a large genomic reference panel. It is based on the Michigan Imputation Server pipeline. Overall, the pipeline filters, phases, and performs imputation on a multi-sample VCF. It outputs the imputed VCF along with key imputation metrics.
