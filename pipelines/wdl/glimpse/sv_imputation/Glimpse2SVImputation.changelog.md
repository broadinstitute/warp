# 0.0.20
2026-08-17 (Date of Last Commit)

* Adds an optional pipeline_header_line input to the WDL

# 0.0.19
2026-08-17 (Date of Last Commit)

* make all task disk sizes dynamic and set bootDiskSizeGb to 0
* update default value of sample_batch_size to 1000

# 0.0.18
2026-08-13 (Date of Last Commit)

* Update `preprocess_pls_gvcf_pipeline_version` to 0.0.10 to pick up removal of background heartbeat monitor logging from merge/concat tasks

# 0.0.17
2026-08-13 (Date of Last Commit)

* Update `preprocess_pls_gvcf_pipeline_version` to 0.0.7 and `batch_pipeline_version` to 0.0.12 for resource optimizations in `PreprocessPLs`, `GLIMPSE2Ligate`, and `PopAndMarginalizeCollisions` tasks

# 0.0.16
2026-08-12 (Date of Last Commit)

* update batch_pipeline_version to 0.0.12 and preprocess_pls_gvcf_pipeline_version to 0.0.8

# 0.0.15
2026-08-12 (Date of Last Commit)

* Replace top-level FOFN inputs with a `gvcf_manifest` input that expects `gvcf_path`, `gvcf_index_path`, and `sample_id` columns to exist.
* Batch samples from the manifest, while still converting array inputs into a manifest internally for compatibility.
* Pass each manifest batch directly into `PreprocessPLsGVCF`, which now parses the manifest internally.
* Use the sample count emitted by `PreprocessPLsGVCF` instead of counting manifest rows in the top-level workflow.

# 0.0.14
2026-08-11 (Date of Last Commit)

* Add optional `info_filter_for_inclusion` input; when set above 0.0, variants with INFO score below the threshold are excluded from the final per-chromosome popped output VCFs.

# 0.0.13
2026-08-09 (Date of Last Commit)

* Add array-input sample batching for `input_gvcfs`, `input_gvcf_idxs`, and `sample_ids`.
* Support FOFN input path (`input_gvcfs_fofn`, `input_gvcf_idxs_fofn`, `sample_ids_file`) as a single batch.
* Run `PreprocessPLsGVCF` and `Glimpse2SVImputationBatch` per sample batch.
* Merge per-chromosome popped outputs across batches and recompute cohort-level `AF` and `INFO`.
* Emit only `glimpse2_popped_posteriors_vcf` outputs from the top-level workflow.

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
