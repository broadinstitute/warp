---
sidebar_position: 7
slug: /All_of_Us/RNA_Seq_QTL/calculate_phenotype_groups
title: Calculate Phenotype Groups
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.0.1](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/CalculatePhenotypeGroups.changelog.md) | January, 2026 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the Calculate Phenotype Groups workflow

[`CalculatePhenotypeGroups.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/CalculatePhenotypeGroups.wdl) generates phenotype grouping metadata for splice-QTL analyses.

It prepares splice data from sample lists and splice matrices and emits a phenotype-groups table used by downstream TensorQTL-style sQTL workflows.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `SampleList` | Sample list for phenotype-group computation. | File |
| `SpliceData` | Splicing data input matrix/file. | File |
| `OutputPrefix` | Prefix used for output naming. | String |
| `memory` | Runtime memory (GB). | Int |
| `disk_space` | Runtime disk size. | Int |
| `num_threads` | Runtime CPU thread count. | Int |

## Outputs

| Output variable name | Filename, if applicable | Output format and description |
| --- | --- | --- |
| `PhenotypeGroups` | `<OutputPrefix>.phenotype_groups.tsv` | Phenotype-group mapping file for sQTL covariate/model workflows. |

## Workflow and WDL

- Workflow: `CalculatePhenotypeGroups`
- Source WDL: [`CalculatePhenotypeGroups.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/CalculatePhenotypeGroups.wdl)

## Versioning

All `CalculatePhenotypeGroups` releases are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/CalculatePhenotypeGroups.changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
