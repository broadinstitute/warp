# 1.0.0
2026-04-15 (Date of Last Commit)

* Initial release of pipeline to perform QC checks on inputs to the Low Pass Imputation pipeline using GLIMPSE2.
* Checks include:
  - If manifest input, all required columns are present
  - Same number of CRAMs, CRAM indices, and sample IDs provided
  - Unique CRAM paths and sample IDs provided
  - Each file input exists and is accessible in an non-requester-pays GCS bucket
  - Each CRAM file size is less than 10GB
