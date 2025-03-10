# 1.0.0
2025-02-26 (Date of Last Commit)

* Initial public release of the ImputationBeagle pipeline.
  * The ImputationBeagle pipeline imputes missing genotypes from a multi-sample VCF using a large genomic reference panel. It is based on the Michigan Imputation Server pipeline but uses the Beagle imputation tool instead of minimac. Overall, the pipeline filters, phases, and performs imputation on a multi-sample VCF. It outputs the imputed VCF.
