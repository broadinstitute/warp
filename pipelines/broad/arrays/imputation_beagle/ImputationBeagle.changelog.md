# 2.0.2
2025-09-03 (Date of Last Commit)

* Add optional pipeline_header_line input that when supplied, will add a header line containing this value to the 
header of the output vcf

# 2.0.1
2025-08-26 (Date of Last Commit)

* Update tasks to use maxRetries of 1 instead of 2. 1 retry is sufficient for transient errors and helps reduce costs.
* Split the two nested scatters from each other so no shard will perform phasing/imputation if any shard fails qc checks
* updated phase/impute tasks to work with the retry with more memory workflow option to help with OOM issues we're seeing in those tasks

# 2.0.0
2025-07-31 (Date of Last Commit)

### Breaking changes
* Update ImputationBeagle pipeline to split by chunks of samples to help scale the workflow to more samples.
This also includes tasks to split the input VCF into sample chunks, run the imputation on each chunk, and then
merge the results back together.  These changes change the scientific output of the pipeline.
* Removed SeparateMultiallelics and RemoveSymbolicAlleles tasks as they are not something we want to do
anymore.

### Additional changes
* Use SSD for tasks that localize large files to help with performance on Google Batch.
* Add maxRetries to all tasks to help with performance on Google Batch.

# 1.0.2
2025-04-07 (Date of Last Commit)

* Update Imputation Tasks to not request external IP addresses (use runtime attribute 'noAddress: true').

# 1.0.1
2025-04-01 (Date of Last Commit)

* Have ImputationBeagle use gatk docker in GAR rather than pull it from dockerhub
* Update Imputation Tasks to use dockers (ubuntu and tidyverse) that have been moved to GAR.

# 1.0.0
2025-02-26 (Date of Last Commit)

* Initial public release of the ImputationBeagle pipeline.
  * The ImputationBeagle pipeline imputes missing genotypes from a multi-sample VCF using a large genomic reference
  * panel. It is based on the Michigan Imputation Server pipeline but uses the Beagle imputation tool instead of minimac. Overall, the pipeline filters, phases, and performs imputation on a multi-sample VCF. It outputs the imputed VCF.
