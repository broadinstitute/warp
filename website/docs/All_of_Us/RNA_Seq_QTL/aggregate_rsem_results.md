---
sidebar_position: 4
slug: /All_of_Us/RNA_Seq_QTL/aggregate_rsem_results
title: RSEM Results Aggregation
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [Changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/GTEx/aggregate_rsem_results.changelog.md) | See changelog | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the RSEM Results Aggregation workflow

[`aggregate_rsem_results.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/GTEx/aggregate_rsem_results.wdl) defines workflow `aggregate_rsem_results`, which combines per-sample RSEM transcript- and gene-level results into cohort-level expression matrices.

The workflow accepts newline-delimited lists of RSEM isoform and gene result paths. It stages each set of files and runs the aggregation script separately to produce transcript-level TPM, isoform percentage, and expected-count matrices, plus gene-level TPM and expected-count matrices.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `rsem_isoforms_list` | Newline-delimited list of per-sample RSEM isoform result file paths. | File |
| `rsem_genes_list` | Newline-delimited list of per-sample RSEM gene result file paths. | File |
| `prefix` | Prefix used for the aggregated output filenames. | String |
| `memory` | Memory allocated to the aggregation task, in GB. | Int |
| `disk_space` | Disk space allocated to the aggregation task, in GB. | Int |
| `num_threads` | Number of CPUs allocated to the aggregation task. | Int |
| `num_preempt` | Number of preemptible retries for the aggregation task. | Int |
| `pipeline_version` | Version string recorded in the workflow output. Defaults to `aou_9.0.1`. | String |

## Outputs

| Output variable name | Filename, if applicable | Output format and description |
| --- | --- | --- |
| `transcripts_tpm` | `<prefix>.rsem_transcripts_tpm.txt.gz` | Cohort-level transcript TPM matrix. |
| `transcripts_isopct` | `<prefix>.rsem_transcripts_isopct.txt.gz` | Cohort-level transcript isoform-percentage matrix. |
| `transcripts_expected_count` | `<prefix>.rsem_transcripts_expected_count.txt.gz` | Cohort-level transcript expected-count matrix. |
| `genes_tpm` | `<prefix>.rsem_genes_tpm.txt.gz` | Cohort-level gene TPM matrix. |
| `genes_expected_count` | `<prefix>.rsem_genes_expected_count.txt.gz` | Cohort-level gene expected-count matrix. |
| `pipeline_version_out` | N/A | Pipeline version string supplied through `pipeline_version`. |

## Workflow and WDL

- Workflow: `aggregate_rsem_results`
- Source WDL: [`aggregate_rsem_results.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/GTEx/aggregate_rsem_results.wdl)

## Processing stages

1. **Aggregate transcript results:** downloads the isoform result files and generates transcript TPM, isoform-percentage, and expected-count matrices.
2. **Aggregate gene results:** downloads the gene result files and generates gene TPM and expected-count matrices.
3. **Return outputs:** provides the five compressed matrices and the pipeline version string.

## Versioning

Releases of `aggregate_rsem_results` are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/GTEx/aggregate_rsem_results.changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>