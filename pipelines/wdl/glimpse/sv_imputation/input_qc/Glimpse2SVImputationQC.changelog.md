# 0.0.1
2026-08-24 (Date of Last Commit)

* Initial release of pipeline to perform QC checks on inputs to the SV Imputation pipeline using GLIMPSE2.
* Checks include:
  - All required manifest columns (`gvcf_path`, `gvcf_index_path`) are present
  - Unique GVCF paths provided
  - GVCF and GVCF index files have the expected extensions and matching basenames
  - Each GVCF and GVCF index file exists and is accessible in a non-requester-pays GCS bucket
  - Each GVCF file size is less than the configured maximum (default 10GB)
  - Each GVCF header contains the expected chromosomes and chromosome lengths, based on the corresponding reference dictionary
  - Each GVCF header contains FORMAT field definitions for PL and GT
  - Each GVCF header contains multiple GVCFBlocks lines (indicating true GVCF format)
  
