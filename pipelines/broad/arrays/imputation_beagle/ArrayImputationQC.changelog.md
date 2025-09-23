# 1.1.0
2025-09-23 (Date of Last Commit)

* Allow subset of chromosomes in VCF input, but require at least one of the specified contigs.

# 1.0.3
2025-09-15 (Date of Last Commit)

* Add instruction for updating ImputationBeagle wdl when this workflow's pipeline_version changes for improved version tracking across wdls

# 1.0.2
2025-09-03 (Date of Last Commit)

* Add optional pipeline_header_line input to match beagle imputation pipeline inputs

# 1.0.1
2025-08-26 (Date of Last Commit)

* Update tasks to use maxRetries of 1 instead of 2. 1 retry is sufficient for transient errors and helps reduce costs.

# 1.0.0
2025-08-11 (Date of Last Commit)

* Initial release of pipeline to perform QC checks on inputs of the Imputation Beagle pipeline.
* Checks include:
  - Called against HG38 reference genome
  - VCF version 4.x 
  - All chromosomes are present
  - Not a whole genome sequencing (WGS) file
