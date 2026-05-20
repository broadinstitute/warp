# 1.0.2
2026-05-20 (Date of Last Commit)

* Add check for correct reference alignment in input CRAMs 

# 1.0.1
2026-05-11 (Date of Last Commit)

* Use cram_manifest if provided instead of crams/cram_indices/sample_ids input arrays

# 1.0.0
2026-04-15 (Date of Last Commit)

* Initial release of pipeline to perform QC checks on inputs to the Low Pass Imputation pipeline using GLIMPSE2.
* Checks include:
  - If manifest input, all required columns are present
  - Same number of CRAMs, CRAM indices, and sample IDs provided
  - Unique CRAM paths and sample IDs provided
  - Each file input exists and is accessible in an non-requester-pays GCS bucket
  - Each CRAM file size is less than 10GB
