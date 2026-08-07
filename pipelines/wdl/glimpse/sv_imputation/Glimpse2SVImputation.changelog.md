# 0.0.13
2026-08-07 (Date of Last Commit)

* Added optional `info_filter_for_inclusion` input to filter output VCFs by INFO score
* Updated batch_pipeline_version to 0.0.11

# 0.0.12
2026-08-06 (Date of Last Commit)

* Update batch_pipeline_version to 0.0.10

# 0.0.11
2026-08-05 (Date of Last Commit)

* update batch_pipeline_version

# 0.0.10
2026-08-04 (Date of Last Commit)

* Update Glimpse2 docker image to tag `imputation-glimpse2:1.2.0-8671138-1784681771`
* Update `batch_pipeline_version` to 0.0.8
* Update `preprocess_pls_gvcf_pipeline_version` to 0.0.6

# 0.0.9
2026-08-03 (Date of Last Commit)

* Updated batch_pipeline_version to 0.0.7

# 0.0.8
2026-07-31 (Date of Last Commit)

* update batch_pipeline_version

# 0.0.7
2026-07-30 (Date of Last Commit)

* rename entity_ids input to sample_ids
* remove sample_names_map_file input as preprocess no longer uses it

# 0.0.6
2026-07-23 (Date of Last Commit)

* Update workflow to process array of chromosomes instead of just one
* remove now unused remapping input

# 0.0.5
2026-07-21 (Date of Last Commit)

* Updated both `preprocess_pls_gvcf_pipeline_version` and `batch_pipeline_version` to 0.0.4

# 0.0.4
2026-07-17 (Date of Last Commit)

* * Updated tasks to use the official bcftools-vcftools and sv-imputation-rust-tools docker images
* sv-imputation-rust-tools contains the 3 rust binaries and tasks have been updated accordingly
  to use the binaries from this image instead of building them in the pipeline

# 0.0.3
2026-07-16 (Date of Last Commit)

* Optional cpu override input for Glimpse Phase

# 0.0.2
2026-07-14 (Date of Last Commit)

* remove output_prefix input from PreProcessGVCFs call

# 0.0.1
2026-07-08 (Date of Last Commit)

* Early draft of the glimpse sv imputation pipeline (putting this here to satisfy PR checks for now)
