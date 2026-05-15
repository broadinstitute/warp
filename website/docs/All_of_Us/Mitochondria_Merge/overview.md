---
title: Mitochondria Merge Pipeline Overview
sidebar_position: 1
slug: /All_of_Us/Mitochondria_Merge/README
className: aou-doc-page
---

<div className="aou-folder-text">

# Mitochondria Merge Pipeline Overview

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.1.0](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/merge/mitochondria_merge.changelog.md) | May, 2026 | [WARP Pipelines](mailto:warp@broadinstitute.org) | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the Mitochondria Merge workflow

The mtDNA Coverage Merge Pipeline ([`mitochondria_merge.wdl`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/merge/mitochondria_merge.wdl)) is a cloud-optimized WDL pipeline that takes per-sample mtDNA VCF files and coverage metrics produced by the [Mitochondria Single Sample](../Mitochondria_Single_Sample/README.md) and merges them into a single annotated cohort-wide callset. The pipeline is designed to scale to hundreds of thousands of samples and has been validated on the All of Us (AoU) v9 dataset of 535,000 samples.

The pipeline proceeds in six stages: metadata preprocessing, HDF5 coverage database construction, sharded VCF ingestion and multi-round merging, sharded homoplasmic-reference imputation and artifact filtering, and finally a dual-track annotation step that produces both a complete callset and a QC-filtered callset.

A companion utility pipeline, [`get_wgs_median_coverage.wdl`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/merge/get_wgs_median_coverage.wdl), can be used to extract WGS median coverage values from per-sample Picard metrics files and produce the `wgs_median_coverage_tsv` input required by this pipeline.

## Quickstart table

| Pipeline Feature | Description |
| :--: | :-- |
| Assay type | Whole-genome sequencing (mtDNA, cohort-level) |
| Overall workflow | Metadata join, coverage DB construction, sharded VCF merge, hom-ref imputation, annotation |
| Workflow language | WDL 1.0 |
| Genomic reference sequence | hg38 |
| Data input file format | TSV (per-sample VCF and coverage file paths) |
| Data output file format | Hail MatrixTable (tar.gz), VCF, TSV |

## Set-up

### Mitochondria Merge installation and requirements

The pipeline code can be downloaded by cloning the [WARP GitHub repository](https://github.com/broadinstitute/warp). The WDL is located at [`all_of_us/mitochondria/mitochondria_merge.wdl`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/merge/mitochondria_merge.wdl).

The pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant workflow management system.

### Preparing inputs with get_wgs_median_coverage

The [`get_wgs_median_coverage.wdl`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/merge/get_wgs_median_coverage.wdl) utility pipeline reads a two-column TSV (`research_id`, `metrics_path`) pointing to per-sample Picard WGS metrics files and extracts the median coverage value for each sample. Its output (`median_coverage_tsv`) is the `wgs_median_coverage_tsv` input to `mitochondria_merge`. The pipeline shards the input TSV for parallel processing and merges results at the end.

## Inputs

### Sample data inputs

| Input variable name | Description | Type |
| --- | --- | --- |
| `full_data_tsv` | Per-sample data table with paths to VCF and coverage files from the Mitochondria Single Sample Pipeline | File |
| `sample_list_tsv` | *(Optional)* Single-column list of sample IDs to subset to | File? |

### Metadata inputs

| Input variable name | Description | Type |
| --- | --- | --- |
| `coverage_tsv` | Per-sample QC metrics (mean WGS coverage, contamination, biosample collection date) | File |
| `ancestry_tsv` | Predicted genetic ancestry per sample | File |
| `dob_tsv` | Date of birth per sample | File |
| `wgs_median_coverage_tsv` | WGS median coverage per sample (produced by `get_wgs_median_coverage.wdl`) | File |

### GCS bucket inputs

| Input variable name | Description | Type |
| --- | --- | --- |
| `vcf_merge_output_bucket` | GCS bucket for intermediate MatrixTable tarballs (VCF shards, merged MTs, finalized shards) | String |
| `annotated_output_bucket` | GCS bucket for final annotated outputs | String |

### Optional tuning inputs

| Input variable name | Default | Description | Type |
| --- | --- | --- | --- |
| `vcf_col_name` | `"final_vcf"` | Column name in `full_data_tsv` containing the per-sample VCF path | String |
| `combined_mt_name` | `"combined_vcf"` | Name for the combined output MatrixTable | String |
| `vcf_merge_shard_size` | `2500` | Samples per VCF ingest shard | Int |
| `vcf_merge_merge_fanin` | `10` | Fan-in per merge round | Int |
| `vcf_merge_shard_n_partitions` | `192` | Hail partitions per shard MT | Int |
| `finalize_shard_size` | `25000` | Samples per finalize shard | Int |
| `finalize_shard_n_partitions` | `256` | Hail partitions per finalize shard MT | Int |
| `finalize_union_n_partitions` | `1000` | Hail partitions for the final unioned MT | Int |

## Mitochondria Merge tasks and tools

The pipeline calls a series of tasks defined in [`mitochondria_merge.wdl`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/merge/mitochondria_merge.wdl). The following sections describe each step:

1. [(Optional) Subset data table](#1-optional-subset-data-table)
2. [Preprocess metadata TSV](#2-preprocess-metadata-tsv)
3. [Build HDF5 coverage database](#3-build-hdf5-coverage-database)
4. [Shard, ingest, and merge VCFs](#4-shard-ingest-and-merge-vcfs)
5. [Finalize: impute hom-ref and filter artifacts](#5-finalize-impute-hom-ref-and-filter-artifacts)
6. [Add annotations](#6-add-annotations)

| Task name | Docker image | Description |
| --- | --- | --- |
| [`subset_data_table`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/merge/mitochondria_merge.wdl) | `python-numpy-pandas:1.0.0-2.2.3-1.25.2` | *(Optional)* Filters the input data table to a specified sample list |
| [`process_tsv_files`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/merge/mitochondria_merge.wdl) | `python-numpy-pandas:1.0.0-2.2.3-1.25.2` | Joins metadata TSVs into a master sample sheet |
| [`annotate_coverage`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/merge/mitochondria_merge.wdl) | `aou-mito-coverage-db:1.0.0` | Builds an HDF5 coverage database from per-sample coverage TSVs |
| [`make_vcf_shards_from_tsv`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/merge/mitochondria_merge.wdl) | `python-numpy-pandas:1.0.0-2.2.3-1.25.2` | Partitions the sample list into shard TSVs |
| [`build_vcf_shard_mt`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/merge/mitochondria_merge.wdl) | `aou-mito-hail-processing:1.0.0` | Imports one shard of VCFs into a Hail MatrixTable |
| [`make_mt_merge_groups`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/merge/mitochondria_merge.wdl) | `python-numpy-pandas:1.0.0-2.2.3-1.25.2` | Plans fan-in merge groups from a list of MT tarballs |
| [`merge_mt_shards`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/merge/mitochondria_merge.wdl) | `aou-mito-hail-processing:1.0.0` | Merges a group of MT tarballs into one MT via `multi_way_union_mts` |
| [`shard_mt_by_samples`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/merge/mitochondria_merge.wdl) | `aou-mito-hail-processing:1.0.0` | Splits the merged MT into per-sample-shard MT tarballs |
| [`finalize_mt_with_covdb`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/merge/mitochondria_merge.wdl) | `aou-mito-hail-processing:1.0.0` | Imputes hom-ref genotypes from HDF5 coverage DB and applies artifact filters |
| [`union_mt_shards`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/merge/mitochondria_merge.wdl) | `aou-mito-hail-processing:1.0.0` | Unions finalized sample-shard MTs into a single cohort MT |
| [`add_annotations`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/merge/mitochondria_merge.wdl) | `aou-mito-hail-processing:1.0.0` | Adds cohort-wide variant statistics, QC filters, and population/haplogroup allele frequencies |

### 1. (Optional) Subset data table

`subset_data_table` filters `full_data_tsv` to only the rows matching sample IDs in `sample_list_tsv`. If no sample list is provided, the full table passes through unchanged.

### 2. Preprocess metadata TSV

`process_tsv_files` joins the sample data table with `coverage_tsv`, `ancestry_tsv`, `dob_tsv`, and `wgs_median_coverage_tsv` on sample ID using pandas. Columns are renamed to the expected downstream format (`ancestry_pred` → `pop`, etc.) and an `age` field is derived from biosample collection date and date of birth. The output `processed_tsv` is the master sample sheet for all downstream tasks.

### 3. Build HDF5 coverage database

`annotate_coverage` reads the per-sample coverage TSVs listed in `processed_tsv` and writes coverage values into an HDF5 file (`coverage.h5`) with shape `(n_samples, n_positions)`. Per-position summary statistics (mean, median, fraction of samples with >100× and >1000× coverage) are computed using numpy. This step runs entirely without Hail or Spark.

Output `coverage_db.tar.gz` contains `coverage.h5` and `coverage_summary.tsv` and is passed to the finalize and annotation tasks.

### 4. Shard, ingest, and merge VCFs

VCF ingestion is parallelized across sample shards to avoid driver memory exhaustion at large cohort sizes:

- **`make_vcf_shards_from_tsv`** partitions the sample list into TSV shards of `vcf_merge_shard_size` samples each.
- **`build_vcf_shard_mt`** (scattered) imports one shard of per-sample VCFs into a Hail MatrixTable and writes the result as a `.tar.gz` to `vcf_merge_output_bucket`.
- **`make_mt_merge_groups` + `merge_mt_shards`** (3 rounds, scattered) merge shard MTs in a fan-in tree with `vcf_merge_merge_fanin` MTs per group. Three rounds reduce hundreds of shard MTs to a single merged MT. Intermediate tarballs are persisted in `vcf_merge_output_bucket` for call-caching and restartability.

### 5. Finalize: impute hom-ref and filter artifacts

The merged MT is further sharded by samples for parallelized finalization:

- **`shard_mt_by_samples`** splits the merged MT into sample shards of `finalize_shard_size` samples each.
- **`finalize_mt_with_covdb`** (scattered) processes each sample shard: reads the HDF5 coverage DB in position blocks, filters to only the rows in each block before annotating entries (avoiding repeated full-matrix scans), imputes homoplasmic-reference genotypes where coverage exceeds the threshold, and applies the artifact-prone site filter.
- **`union_mt_shards`** reunites the finalized sample shards into a single cohort-wide MT via `union_cols`.

### 6. Add annotations

`add_annotations` runs `add_annotations.py` on the final unioned MT to compute cohort-wide variant statistics, per-population and per-haplogroup allele frequencies, and QC filters. It is called twice in parallel:

- `annotated` (`keep_all_samples = true`): all samples retained
- `filt_annotated` (`keep_all_samples = false`): samples failing QC filters (contamination, copy number, etc.) removed

Outputs are written directly to timestamped subdirectories under `annotated_output_bucket`.

## Outputs

| Output variable name | Description | Type |
| --- | --- | --- |
| `processed_tsv` | Master sample sheet after metadata joins | File |
| `output_coverage_db` | `coverage_db.tar.gz` — HDF5 coverage database and summary TSV | File |
| `combined_vcf` | Tar.gz of the finalized, unioned cohort MatrixTable | File |
| `annotated_output_gcs_path` | GCS path to annotated outputs (all samples) | String |
| `filt_annotated_output_gcs_path` | GCS path to annotated outputs (QC-filtered samples only) | String |

## Downstream pipeline: Mito Post Processing

After `MitochondriaMerge` completes, the annotated Hail MatrixTable can be passed to the optional [`mito_post_processing.wdl`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/merge/mito_post_processing.wdl) pipeline, which filters the callset and generates a standard set of QC plots.

### Mito Post Processing inputs

| Input variable name | Description | Type |
| --- | --- | --- |
| `input_path` | GCS path to the annotated Hail MatrixTable (output of `add_annotations`) | String |
| `output_path` | GCS destination directory for all outputs | String |
| `output_base` | Filename prefix applied to every output file | String |
| `hail_docker` | Docker image (default: `us.gcr.io/broad-gotc-prod/aou_mitochondria_post:0.0.5`) | String |

### Mito Post Processing task

The pipeline has a single task, `RunMitoPostProcessing`, which runs `mito_plot_filter.py` to apply variant and sample filters and produce all outputs.

| Task name | Docker image | Description |
| --- | --- | --- |
| [`RunMitoPostProcessing`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/merge/mito_post_processing.wdl) | `aou_mitochondria_post:0.0.5` | Filters the annotated MatrixTable and generates QC plots |

### Mito Post Processing outputs

| Output variable name | Description | Type |
| --- | --- | --- |
| `filtered_vcf` | bgzipped VCF of the filtered callset | File |
| `filtered_vcf_tbi` | Tabix index for `filtered_vcf` | File |
| `sample_metadata_tsv` | Per-sample metadata table | File |
| `variants_per_sample_svg` | Distribution of variant counts per sample | File |
| `mito_cn_distribution_svg` | mtDNA copy number distribution across samples | File |
| `variant_allele_frequency_svg` | Variant allele frequency spectrum | File |
| `variant_af_and_allele_fraction_svg` | Combined AF and allele fraction plot | File |
| `numt_fp_by_mtcn_svg` | NuMT false-positive rate as a function of mtCN | File |
| `haplogroup_heteroplasmy_svg` | Heteroplasmy distribution per haplogroup | File |
| `haplogroup_homoplasmy_svg` | Homoplasmy distribution per haplogroup | File |

## Versioning and testing

All mitochondria merge pipeline releases are documented in the [mitochondria merge changelog](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/merge/mitochondria_merge.changelog.md).

## Citing the Mitochondria Merge Pipeline

When citing this pipeline, please cite the original mtSwirl publication:

Gupta, R., Kanai, M., Durham, T.J. et al. Nuclear genetic control of mtDNA copy number and heteroplasmy in humans. *Nature*, 2023. https://doi.org/10.1038/s41586-023-06426-5.

When citing WARP, please use:

Kylee Degatano, Aseel Awdeh, Robert Sidney Cox III, Wes Dingman, George Grant, Farzaneh Khajouei, Elizabeth Kiernan, Kishori Konwar, Kaylee L Mathews, Kevin Palis, Nikelle Petrillo, Geraldine Van der Auwera, Chengchen (Rex) Wang, Jessica Way. "Warp Analysis Research Pipelines: Cloud-optimized workflows for biological data processing and reproducible analysis." *Bioinformatics*, 2025; https://doi.org/10.1093/bioinformatics/btaf494.

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues).
