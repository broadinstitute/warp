---
sidebar_position: 7
slug: /All_of_Us/Ancestry_Analysis/run_sample_outlier_qc_plotting
title: Run Sample Outlier QC Plotting
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.0.0](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_sample_outlier_qc_plotting.changelog.md) | September, 2025 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the Run Sample Outlier QC Plotting workflow

[`run_sample_outlier_qc_plotting`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_sample_outlier_qc_plotting.wdl) is a WDL workflow that joins sample-level QC/ancestry annotations with demographics and generates interactive visualization outputs.

The workflow reads Hail table artifacts from `run_sample_outlier_qc`, enriches records with demographic race/ethnicity labels, then produces two HTML reports: a principal component scatter plot and a multi-tab QC metric/fitting visualization.

## Quickstart table

| Pipeline Feature | Description | Source |
| :--: | :-- | :--: |
| Analysis type | Demographic join + ancestry/QC visualization | |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Data input file format | Demographics TSV + tarred Hail table from outlier QC | |
| Data output file format | Interactive HTML plots | |
| Primary software | Hail + bokeh | [Hail](https://hail.is/), [Bokeh](https://bokeh.org/) |

## Set-up

### Run Sample Outlier QC Plotting installation and requirements

The workflow code can be downloaded by cloning the [WARP GitHub repository](https://github.com/broadinstitute/warp). For the latest release, please see the [run_sample_outlier_qc_plotting changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_sample_outlier_qc_plotting.changelog.md).

The pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant workflow management system.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `aou_demographics_tsv` | Demographics TSV containing `research_id` and race/ethnicity metadata fields. | File |
| `output_prefix` | Prefix applied to visualization outputs. | String |
| `ancestry_with_flagged_samples_tar_gz` | Tarred Hail table output from `run_sample_outlier_qc` (`<input_prefix>.full.ht.tar.gz`). | File |
| `input_prefix` | Prefix used for input Hail table names generated upstream. | String |

## Run Sample Outlier QC Plotting tasks and tools

The workflow joins demographics and produces two types of interactive reports.

1. [Join demographics to ancestry/QC table](#1-join-demographics-to-ancestryqc-table)
2. [Generate PC scatter plot](#2-generate-pc-scatter-plot)
3. [Generate QC metric fitting report](#3-generate-qc-metric-fitting-report)

To see specific tool parameters, select the task WDL link in the table; then view the `command {}` section of the task in the WDL script.

| Task name and WDL link | Tool | Software | Description |
| --- | --- | --- | --- |
| [join_ancestry_to_demographics](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_sample_outlier_qc_plotting.wdl) | Hail table join | `hailgenetics/hail:0.2.67` | Joins demographics data into ancestry/QC table and writes tarred Hail artifact. |
| [plot_first_pcs](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_sample_outlier_qc_plotting.wdl) | Hail plotting | `hailgenetics/hail:0.2.67` | Generates interactive PC1-vs-PC2 ancestry plot. |
| [plot_metrics_and_fitting](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_sample_outlier_qc_plotting.wdl) | Hail + bokeh tabs | `hailgenetics/hail:0.2.67` | Generates interactive multi-metric QC fitting visualization. |

### 1. Join demographics to ancestry/QC table

`join_ancestry_to_demographics` extracts the upstream Hail table, joins records by sample ID to demographics metadata, and writes `<output_prefix>.ancestry_with_flagged_samples_demographics.ht.tar.gz`.

### 2. Generate PC scatter plot

`plot_first_pcs` creates `<output_prefix>.pc1vspc2.html`, visualizing ancestry clusters and overlaying self-reported race/ethnicity metadata in hover fields.

### 3. Generate QC metric fitting report

`plot_metrics_and_fitting` creates `<output_prefix>.metrics.html`, an interactive tabbed report showing ancestry-stratified metric distributions and fitted trends.

## Outputs

| Output variable name | Filename, if applicable | Output format and description |
| --- | --- | --- |
| `pcs_plot` | `<output_prefix>.pc1vspc2.html` | Interactive PC1-vs-PC2 ancestry scatter plot with demographic hover metadata. |
| `metrics_plot` | `<output_prefix>.metrics.html` | Interactive tabbed QC metric and fitting visualization across selected QC metrics. |

## Versioning

All `run_sample_outlier_qc_plotting` releases are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_sample_outlier_qc_plotting.changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
