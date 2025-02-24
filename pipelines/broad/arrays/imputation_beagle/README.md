### ImputationBeagle summary

The ImputationBeagle pipeline imputes missing genotypes from a multi-sample VCF using the [Beagle imputation tool](https://faculty.washington.edu/browning/beagle/beagle.html) and a large genomic reference panel. Overall, the pipeline filters, phases, and performs imputation on a multi-sample VCF. 

### ArrayImputationQuotaConsumed summary

The ArrayImputationQuotaConsumed pipeline is used by the All of Us/AnVIL Imputation Service and calculates the number of samples in the input multi-sample VCF, which is the metric used by the service for ImputationBeagle pipeline quota.
