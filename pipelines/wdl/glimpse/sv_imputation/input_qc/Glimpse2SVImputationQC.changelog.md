# 0.0.1
2026-08-24 (Date of Last Commit)

* Initial release of pipeline to perform QC checks on inputs to the SV Imputation pipeline using GLIMPSE2.
* Checks include:
  - If manifest input, all required columns are present
  - Unique GVCF paths and sample IDs provided
  - GVCF and GVCF index files have the expected extensions and matching basenames
  - Each GVCF and GVCF index file exists and is accessible in a non-requester-pays GCS bucket
  - Each GVCF file size is less than the configured maximum
  - Each GVCF header contains the expected chromosomes, with reference alignment MD5 checked against the provided reference dictionary when available
