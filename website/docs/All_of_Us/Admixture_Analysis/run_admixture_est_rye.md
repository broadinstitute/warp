---
sidebar_position: 3
slug: /All_of_Us/Admixture_Analysis/run_admixture_est_rye
title: Admixture Rye
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.0.0](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/run_admixture_est_rye.changelog.md) | September, 2025 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the run_admixture_est_rye workflow

[`run_admixture_est_rye`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/run_admixture_est_rye.wdl) is a WDL workflow that runs the [Rye admixture estimation tool](https://github.com/healthdisparities/rye) to estimate global ancestry proportions from PCA-derived data. It is the second step in the **Admixture Rye** analysis path used in All of Us processing.

The workflow accepts Rye-format eigenvalues, eigenvectors, and a population-to-group mapping file produced by [`run_preprocess_admixture_est_rye`](./run_preprocess_admixture_est_rye.md) and runs the Rye R tool to produce ancestry proportion estimates. Outputs are analogous in format to the `.Q` and `.fam` files produced by the classic ADMIXTURE tool, enabling direct comparison between methods.

Rye uses an iterative optimization approach parameterized by the number of principal components, optimization rounds, and iterations, balancing accuracy against runtime.

## Quickstart table

| Workflow Feature | Description | Source |
| :--: | :-- | :--: |
| Analysis type | Rye-based admixture estimation | |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Input data type | Rye-format eigenvalues, eigenvectors, pop2group mapping | |
| Output file format | `.Q` (ancestry proportions), `.fam` (sample metadata) | |
| Primary tool | Rye | [Rye GitHub](https://github.com/healthdisparities/rye) |
| Part of analysis path | Admixture Rye (Step 2 of 2) | |

## Set-up

### run_admixture_est_rye installation and requirements

The workflow code can be downloaded by cloning the [WARP GitHub repository](https://github.com/broadinstitute/warp). For the latest release, please see the [run_admixture_est_rye changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/run_admixture_est_rye.changelog.md).

The pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant workflow management system.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `eigenvalues_file` | Rye-compatible eigenvalues file produced by [`run_preprocess_admixture_est_rye`](./run_preprocess_admixture_est_rye.md). | File |
| `eigenvec_file` | Rye-compatible eigenvectors file containing merged training and testing PCA projections. | File |
| `pop2group_file` | Population-to-group mapping file used by Rye to define reference ancestry groups. | File |
| `prefix` | Output file prefix (e.g., `aou_delta`). Prepended to all output filenames. | String |
| `pcs` | Number of principal components to use in estimation. Default: `20`. | Int |
| `rounds` | Number of optimization rounds. Higher values increase accuracy but increase runtime. Default: `200`. | Int |
| `iter` | Number of iterations per optimization round. Default: `100`. | Int |
| `cpus` | Number of CPU threads. Default: `16`. | Int |
| `docker_image` | Docker image with Rye installed. Default: `us-central1-docker.pkg.dev/broad-dsde-methods/aou-auxiliary/rye-admixture-estimation-tool:v1.0`. | String |

## run_admixture_est_rye tasks and tools

This workflow calls a single task to run Rye ancestry estimation.

1. [Run Rye admixture estimation](#1-run-rye-admixture-estimation)

To see specific tool parameters, select the task WDL link in the table; then view the `command {}` section of the task in the WDL script.

| Task name and WDL link | Tool | Software | Description |
| --- | --- | --- | --- |
| [run_rye](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/run_admixture_est_rye.wdl) | Rye | R | Runs the Rye R script to estimate ancestry proportions from PCA data. Outputs `.Q` and `.fam` files with admixture proportions and sample metadata. |

### 1. Run Rye admixture estimation

The `run_rye` task invokes the Rye R script (`/rye/rye.R`) with the eigenvalues, eigenvectors, and population mapping inputs. Rye performs iterative optimization to estimate the proportion of ancestry each sample derives from each reference population group. After the Rye run completes, the output `.Q` file is renamed to follow the pattern `<prefix>-<pcs>.Q` for consistent delocalization. See the [Rye documentation](https://github.com/healthdisparities/rye) for details on file formats and examples.

## Outputs

| Output variable name | Filename | Output format and description |
| --- | --- | --- |
| `qFile` | `<prefix>-<pcs>.Q` | Tab-delimited admixture proportion matrix. Each row is a sample; each column is an ancestry component. Analogous to ADMIXTURE `.Q` output. |
| `famFile` | `<prefix>-<pcs>.fam` | Sample metadata file in PLINK `.fam` format. |

## Versioning

All `run_admixture_est_rye` releases are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/run_admixture_est_rye.changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
