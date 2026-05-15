---
sidebar_position: 6
slug: /All_of_Us/Ancestry_Analysis/run_sample_outlier_qc
title: Run Sample Outlier QC
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.0.0](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_sample_outlier_qc.changelog.md) | September, 2025 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the Run Sample Outlier QC workflow

[`run_sample_outlier_qc`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_sample_outlier_qc.wdl) is a WDL workflow that joins ancestry predictions with cohort callset summary statistics and identifies outlier samples using ancestry-stratified quality-control residuals.

The workflow first merges ancestry outputs with requested QC metrics, then computes residual-based filters using gnomAD QC utilities. It outputs both full-sample QC annotations and a filtered table containing only flagged samples.

## Quickstart table

| Pipeline Feature | Description | Source |
| :--: | :-- | :--: |
| Analysis type | Sample-level outlier QC with ancestry-aware stratification | |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Data input file format | TSV/CSV ancestry and callset summary inputs | |
| Data output file format | TSV + tarred Hail tables (`.ht.tar.gz`) | |
| Primary software | Hail + gnomAD QC utilities | [Hail](https://hail.is/), [gnomAD methods](https://github.com/broadinstitute/gnomad_methods) |

## Set-up

### Run Sample Outlier QC installation and requirements

The workflow code can be downloaded by cloning the [WARP GitHub repository](https://github.com/broadinstitute/warp). For the latest release, please see the [run_sample_outlier_qc changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_sample_outlier_qc.changelog.md).

The pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant workflow management system.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `callset_summary_csv` | Cohort callset summary statistics CSV (e.g., GVS summary metrics). | File |
| `ancestry_results_tsv` | Ancestry prediction table output from `run_ancestry.wdl`. | File |
| `output_prefix` | Prefix applied to all output artifacts. | String |
| `metrics_to_check_in` | *(Optional)* Python-list string of metrics to evaluate. Defaults to a preset list of variant/QC metrics. | String? |

## Run Sample Outlier QC tasks and tools

The workflow runs two tasks to create a full ancestry+QC table and then identify outliers.

1. [Join ancestry and summary metrics](#1-join-ancestry-and-summary-metrics)
2. [Compute stratified outlier filters](#2-compute-stratified-outlier-filters)

To see specific tool parameters, select the task WDL link in the table; then view the `command {}` section of the task in the WDL script.

| Task name and WDL link | Tool | Software | Description |
| --- | --- | --- | --- |
| [join_ancestry_to_stats](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_sample_outlier_qc.wdl) | Hail table joins | `hailgenetics/hail:0.2.67` | Joins ancestry results with callset summary metrics and writes full ancestry Hail table artifact. |
| [determine_outlier_qc](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_sample_outlier_qc.wdl) | gnomAD residual filtering | `hailgenetics/hail:0.2.67` | Computes residuals/threshold filters and exports flagged and full-sample outputs. |

### 1. Join ancestry and summary metrics

`join_ancestry_to_stats` imports ancestry and callset summary tables, joins selected metrics by sample ID, and writes `<output_prefix>.full_ancestry.ht.tar.gz`.

### 2. Compute stratified outlier filters

`determine_outlier_qc` computes ancestry-stratified residual metrics, applies QC thresholds, and exports full and flagged sample tables as TSV + Hail artifacts.

## Outputs

| Output variable name | Filename, if applicable | Output format and description |
| --- | --- | --- |
| `flagged_samples_tsv` | `<output_prefix>.flagged_samples.tsv` | TSV of samples with one or more QC metric filter failures. |
| `all_samples_tsv` | `<output_prefix>.full.tsv` | TSV containing full sample set with QC residual/filter annotations. |
| `ancestry_with_flagged_samples_tar_gz` | `<output_prefix>.full.ht.tar.gz` | Tarred Hail table containing full ancestry/QC annotations for all samples. |

## Versioning

All `run_sample_outlier_qc` releases are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_sample_outlier_qc.changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
