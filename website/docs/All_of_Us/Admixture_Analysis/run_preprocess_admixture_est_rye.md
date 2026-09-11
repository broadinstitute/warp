---
sidebar_position: 2
slug: /All_of_Us/Admixture_Analysis/run_preprocess_admixture_est_rye
title: Admixture Rye Preprocessing
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.0.0](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/run_preprocess_admixture_est_rye.changelog.md) | September, 2025 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the run_preprocess_admixture_est_rye workflow

[`run_preprocess_admixture_est_rye`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/run_preprocess_admixture_est_rye.wdl) is a WDL workflow that generates the input files required by the [Rye admixture estimation tool](https://github.com/healthdisparities/rye). It is the first step in the **Admixture Rye** analysis path used in All of Us processing.

The workflow consumes PCA outputs produced upstream by the [ancestry pipeline](https://broadinstitute.github.io/warp/docs/All_of_Us/Ancestry_Analysis/overview) — specifically eigenvalues, training PCA projections, and ancestry prediction data — and reshapes them into the file formats expected by Rye. It uses [Hail](https://hail.is/) and [pandas](https://pandas.pydata.org/) to transform and merge training and testing PCA projections and to export a population-to-group mapping file.

This workflow must be run before [`run_admixture_est_rye`](./run_admixture_est_rye.md).

## Quickstart table

| Workflow Feature | Description | Source |
| :--: | :-- | :--: |
| Analysis type | Admixture preprocessing (PCA → Rye inputs) | |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Input data type | PCA eigenvalues, training PCA table, ancestry prediction TSV | |
| Output file format | TSV (eigenvalues, eigenvectors, pop2group) | |
| Primary tool | Hail, pandas | [Hail](https://hail.is/), [pandas](https://pandas.pydata.org/) |
| Part of analysis path | Admixture Rye (Step 1 of 2) | |

## Set-up

### run_preprocess_admixture_est_rye installation and requirements

The workflow code can be downloaded by cloning the [WARP GitHub repository](https://github.com/broadinstitute/warp). For the latest release, please see the [run_preprocess_admixture_est_rye changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/run_preprocess_admixture_est_rye.changelog.md).

The pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant workflow management system.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `eigenvalues_url` | Path to the PCA eigenvalues file (Hail Table format). Computed using [Hail HWE-normalized PCA](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.hwe_normalized_pca). | String |
| `training_pca_url` | Path to the training PCA projections TSV. Computed from reference panel data (e.g., HGDP + 1KG) using Hail PCA. | String |
| `ancestry_data_url` | Path to the All of Us [ancestry prediction file](https://support.researchallofus.org/hc/en-us/articles/4614687617556) containing projected PCA features for test samples. | String |
| `prefix` | Output file prefix (e.g., `aou_delta`). Prepended to all output filenames. | String |
| `cpus` | Number of CPU threads to allocate to the VM. Default: `16`. | Int |
| `docker_image` | Docker image with Hail installed. Default: `hailgenetics/hail:0.2.67`. | String |

## run_preprocess_admixture_est_rye tasks and tools

This workflow calls a single task to transform upstream PCA outputs into Rye-compatible input files.

1. [Preprocess PCA data for Rye](#1-preprocess-pca-data-for-rye)

To see specific tool parameters, select the task WDL link in the table; then view the `command {}` section of the task in the WDL script.

| Task name and WDL link | Tool | Software | Description |
| --- | --- | --- | --- |
| [run_preprocess](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/run_preprocess_admixture_est_rye.wdl) | Hail, pandas | Python | Merges training and testing PCA projections, exports Rye-format eigenvalues, eigenvectors, and population-to-group mapping files. |

### 1. Preprocess PCA data for Rye

The `run_preprocess` task runs a Python script inside the Hail Docker container. It reads eigenvalues, training PCA, and ancestry prediction data; transforms both training and testing PCA projections into a standardized eigenvector format; concatenates them; and writes three output files consumed by [run_admixture_est_rye](./run_admixture_est_rye.md). The population-to-group mapping is derived from distinct population labels in the training data, excluding the `oth` (other) category.

## Outputs

| Output variable name | Filename | Output format and description |
| --- | --- | --- |
| `rye_eigenval` | `<prefix>_rye.eigenvalues` | Tab-separated eigenvalues file in Rye format. |
| `rye_eigenvec` | `<prefix>_rye.eigenvec` | Tab-separated eigenvectors file containing merged training and testing PCA projections in Rye format. |
| `rye_pop2group` | `<prefix>_rye.pop2group` | Tab-separated population-to-group mapping file derived from reference panel labels. |

## Versioning

All `run_preprocess_admixture_est_rye` releases are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/run_preprocess_admixture_est_rye.changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
