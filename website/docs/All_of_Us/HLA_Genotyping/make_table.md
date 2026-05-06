---
sidebar_position: 3
slug: /All_of_Us/HLA_Genotyping/make_table
title: HLA Make Table
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.0.0](https://github.com/broadinstitute/warp/blob/develop/all_of_us/hla/MakeTable_changelog.md) | August, 2025 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the HLA Make Table workflow

[`MakeTable`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/hla/MakeTable.wdl) aggregates per-sample HLA consensus call files into a cohort-level summary table.

It verifies header consistency across samples and emits a unified matrix where each row corresponds to a sample and each gene is represented by two allele columns.

## Quickstart table

| Pipeline Feature | Description | Source |
| :--: | :-- | :--: |
| Analysis type | Cohort-level HLA consensus aggregation | |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Data input file format | Array of per-sample consensus text files + sample IDs | |
| Data output file format | Tab-delimited summary table | |
| Primary software | Bash text processing | |

## Set-up

### HLA Make Table installation and requirements

The workflow code can be downloaded by cloning the [WARP GitHub repository](https://github.com/broadinstitute/warp). For the latest release, please see the [MakeTable changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/hla/MakeTable_changelog.md).

The pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant workflow management system.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `consensus_calls` | Ordered array of per-sample HLA consensus files. | Array[File] |
| `sample_ids` | Ordered array of sample IDs aligned to `consensus_calls`. | Array[String] |

## HLA Make Table tasks and tools

This workflow runs a single task to combine sample files.

1. [Combine per-sample consensus files](#1-combine-per-sample-consensus-files)

| Task name and WDL link | Tool | Software | Description |
| --- | --- | --- | --- |
| [Combine](https://github.com/broadinstitute/warp/blob/develop/all_of_us/hla/MakeTable.wdl) | shell text processing | `continuumio/anaconda:latest` | Builds a shared header, verifies consistency, and writes cohort-wide result table. |

### 1. Combine per-sample consensus files

`Combine` iterates through consensus files, builds/compares gene headers, appends genotype rows, and prepends sample IDs to create the final output matrix.

## Outputs

| Output variable name | Filename, if applicable | Output format and description |
| --- | --- | --- |
| `result` | `result.txt` | Cohort-level HLA summary table with sample IDs and diploid allele columns per gene. |

## Versioning

All `MakeTable` releases are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/hla/MakeTable_changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
