# 0.0.1
2026-03-20 (Date of Last Commit)

* Initial release of pipeline to perform QC checks on inputs to the Low Pass Imputation pipeline using GLIMPSE2.
* Checks include:
  - If manifest input, all required columns are present
  - Same number of CRAMs, CRAIs, and sample IDs provided
  - Unique CRAM paths and sample IDs provided
  - Each CRAM file size is less than 10GB
