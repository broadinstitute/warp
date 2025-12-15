# 2.4.1
2025-12-05 (Date of Last Commit) 

* Update input_qc_version to 1.2.4 to match latest changes in InputQC wdl.

# 2.4.0
2025-11-25 (Date of Last Commit)

* Update Phase and Impute tasks to use a docker image that contains Beagle JAR 
generated from GitHub repo [tmp-sharing/imp-server](https://github.com/tmp-sharing/imp-server/tree/master).

# 2.3.0
2025-11-12 (Date of Last Commit)

* Perform chunk QC on chunks without padding, rather than chunks with padding
* Add a new contigs_info output containing contig-wide metrics. New/changed tasks to support this:
  * Add ExtractUniqueVariantIds task to create a list and count of unique variants from a vcf
  * Change ReferencePanelContig struct to remove bed file and add unique_variants file for reference panel contigs
  * Rename CountVariantsInChunks to CountUniqueVariantIdsInOverlap and change to use unique_variants files rather than input vcfs and reference panel bed files
  * Add CountValidContigChunks task to perform a simpler evaluation of chunk validity
  * Rename StoreChunksInfo task to StoreMetricsInfo and update to include creating a contig-based metrics file in addition to the existing chunk-based metrics file, reporting additional metrics

# 2.2.4
2025-10-17 (Date of Last Commit)

* Update input_qc_version to 1.2.3 to match latest changes in InputQC wdl

# 2.2.3
2025-10-07 (Date of Last Commit)

* Update input_qc_version to 1.2.2 to match latest changes in InputQC wdl

# 2.2.2
2025-10-03 (Date of Last Commit)

* Resolved merge conflicts and reorganize WDL pipelines into unified directory

# 2.2.1
2025-10-01 (Date of Last Commit)

* update input_qc_version to 1.2.1 to match latest changes in InputQC wdl

# 2.2.0
2025-09-29 (Date of Last Commit)

* Add optional min_dr2_for_inclusion input and FilterVcfByDR2 optional task to allow users to specify a minimum DR2 threshold for including imputed variants in the final output VCF. Variants with DR2 below this threshold will be excluded from the final output VCF.

# 2.1.0
2025-09-26 (Date of Last Commit)

* Add CalculateContigsToProcess task to infer chromosomes present in input VCF and filter them against the allowed contigs specified by the repurposed `contigs` input.

# 2.0.4
2025-09-18 (Date of Last Commit)

* Update CountSamples task to use SSD for more efficient localization of large files on Google Batch.
* Update Phase and Impute tasks to use a docker image running on Java 17 instead of Java 8 for better memory management.
* Update LocalizeAndSubsetVcfToRegion task to decompress/compress to get rid of possible empty blocks from the Java 17 update https://broadworkbench.atlassian.net/browse/TSPS-612

# 2.0.3
2025-09-11 (Date of Last Commit)

* Add pipeline_version values for InputQC and QuotaConsumed wdls for improved version tracking across wdls

# 2.0.2
2025-09-03 (Date of Last Commit)

* Add optional pipeline_header_line input that when supplied, will add a header line containing this value to the header of the output vcf

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
