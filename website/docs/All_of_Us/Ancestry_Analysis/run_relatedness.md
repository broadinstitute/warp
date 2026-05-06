---
sidebar_position: 5
slug: /All_of_Us/Ancestry_Analysis/run_relatedness
title: Run Relatedness
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.1.0](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_relatedness.changelog.md) | October, 2025 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the Run Relatedness workflow

[`run_relatedness`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_relatedness.wdl) is a WDL workflow that computes relatedness across cohort samples and identifies a maximal independent set of samples that can be removed to reduce relatedness confounding in downstream analyses.

The workflow submits a Hail-based relatedness job to a Dataproc cluster using a user-provided submission script and PCA-informed inputs, then runs a second task to compute a maximal independent set from related sample pairs. It returns both pairwise relatedness results and a list of flagged samples.

## Quickstart table

| Pipeline Feature | Description | Source |
| :--: | :-- | :--: |
| Analysis type | Cohort relatedness estimation and sample deconfounding | |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Data input file format | VCF + PCA scores + Dataproc submission script | |
| Data output file format | TSV relatedness matrix and flagged sample list | |
| Primary software | Hail + Google Dataproc | [Hail](https://hail.is/), [Dataproc](https://cloud.google.com/dataproc) |

## Set-up

### Run Relatedness installation and requirements

The workflow code can be downloaded by cloning the [WARP GitHub repository](https://github.com/broadinstitute/warp). For the latest release, please see the [run_relatedness changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_relatedness.changelog.md).

The pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant workflow management system.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `vcf_url` | Path to cohort VCF with genotypes used for relatedness. | String |
| `pca_scores_url` | Path to PCA scores corresponding to samples in `vcf_url`. | String |
| `task_identifier` | Identifier used as output filename prefix. | String |
| `statistics` | Statistic to compute. Default: `'kin'`. | String |
| `min_individual_maf` | Minimum individual-specific minor allele frequency. Default: `0.01`. | Float |
| `block_size` | Block matrix size used by the algorithm. Default: `2048`. | Int |
| `min_kinship` | Minimum kinship threshold for reported sample pairs. Default: `0.1`. | Float |
| `min_partitions` | Minimum number of partitions used for optimization. Default: `1200`. | Int |
| `gcs_output_url` | GCS path for relatedness pipeline outputs. | String |
| `executor_cores` | Spark executor core count. | String |
| `driver_cores` | Spark driver core count. | String |
| `executor_memory` | Spark executor memory setting. | String |
| `driver_memory` | Spark driver memory setting. | String |
| `reference_genome` | Reference genome identifier (e.g., `hg38`). | String |
| `max_idle` | Dataproc cluster max idle time in minutes. Default: `60`. | Int |
| `max_age` | Dataproc cluster max age in minutes. Default: `1440`. | Int |
| `num_workers` | Number of Hail Dataproc workers. | Int |
| `gcs_project` | Google Cloud project ID used for Dataproc. | String |
| `gcs_subnetwork_name` | Subnetwork name for Dataproc networking. Default: `'subnetwork'`. | String |
| `submission_script` | Python script submitted to Dataproc to compute relatedness. | File |
| `region` | Dataproc region. Default: `us-central1`. | String |
| `hail_docker` | Docker image for Dataproc orchestration task. | String |
| `hail_docker_maximal_independent_set` | Docker image for maximal independent set task. | String |

## Run Relatedness tasks and tools

The workflow runs two tasks: one for cluster-based relatedness computation and one for maximal independent set filtering.

1. [Compute pairwise relatedness on Dataproc](#1-compute-pairwise-relatedness-on-dataproc)
2. [Compute maximal independent set](#2-compute-maximal-independent-set)

To see specific tool parameters, select the task WDL link in the table; then view the `command {}` section of the task in the WDL script.

| Task name and WDL link | Tool | Software | Description |
| --- | --- | --- | --- |
| [run_relatedness_task](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_relatedness.wdl) | Hail + Dataproc | `us.gcr.io/broad-dsde-methods/lichtens/hail_dataproc_wdl:1.1` | Creates Dataproc cluster, submits relatedness computation job, and copies back relatedness TSV. |
| [run_maximal_independent_set](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_relatedness.wdl) | Hail graph methods | `hailgenetics/hail:0.2.67` | Computes maximal independent set from related pairs and exports flagged samples. |

### 1. Compute pairwise relatedness on Dataproc

`run_relatedness_task` provisions a temporary Dataproc cluster, submits the user-specified script with cohort inputs, and copies `<task_identifier>_relatedness.tsv` from cluster staging storage.

### 2. Compute maximal independent set

`run_maximal_independent_set` reads the relatedness table and applies `hl.maximal_independent_set` to identify a sample subset for removal, exporting the result to `<task_identifier>.relatedness_flagged_samples.tsv`.

## Outputs

| Output variable name | Filename, if applicable | Output format and description |
| --- | --- | --- |
| `relatedness` | `<task_identifier>_relatedness.tsv` | Pairwise relatedness output table containing sample-pair relatedness values. |
| `relatedness_flagged_samples` | `<task_identifier>.relatedness_flagged_samples.tsv` | Table of samples flagged by maximal independent set for potential exclusion. |

## Versioning

All `run_relatedness` releases are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_relatedness.changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
