---
sidebar_position: 8
slug: /All_of_Us/RNA_Seq_QTL/rnaseqc2_aggregate_batched
title: RNA-SeQC2 Batched Aggregation
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.1.0](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/GTEx/rnaseqc2_aggregate_batched.changelog.md) | September, 2026 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the RNA-SeQC2 Batched Aggregation workflow

[`rnaseqc2_aggregate_batched.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/GTEx/rnaseqc2_aggregate_batched.wdl) defines workflow `rnaseqc2_aggregate_batched_workflow`, which aggregates per-sample RNA-SeQC2 outputs into cohort-level files.

The workflow validates a sample manifest, divides the samples into bounded batches, downloads and merges each batch, and performs a final cohort-level merge. This bounded approach limits the number of files processed by an individual task while supporting large cohorts.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `sample_manifest` | Manifest containing the RNA-SeQC2 output paths and sample identifiers to aggregate. | File |
| `prefix` | Prefix used for intermediate and final output filenames. | String |
| `merge_exons` | Whether to merge exon-count GCT files. | Boolean |
| `batch_size` | Maximum number of samples processed in each batch (default `100`). | Int |
| `docker_image` | Container image used by the validation, batch aggregation, and final merge tasks. | String |
| `validation_memory_gb` | Memory allocated to manifest validation in GB (default `1`). | Int |
| `validation_disk_space_gb` | Disk allocated to manifest validation in GB (default `10`). | Int |
| `batch_memory_gb` | Memory allocated to each batch aggregation task in GB (default `4`). | Int |
| `batch_disk_space_gb` | Disk allocated to each batch aggregation task in GB (default `100`). | Int |
| `merge_memory_gb` | Memory allocated to the final merge task in GB (default `8`). | Int |
| `merge_disk_space_gb` | Disk allocated to the final merge task in GB. | Int |
| `num_threads` | Number of threads used for parallel file transfers in batch aggregation (default `8`). | Int |
| `num_preempt` | Number of preemptible retries for each task (default `2`). | Int |

## Outputs

| Output variable name | Filename, if applicable | Output format and description |
| --- | --- | --- |
| `metrics` | `<prefix>.metrics.txt.gz` | Cohort-level aggregated RNA-SeQC2 metrics. |
| `tpm_gct` | `<prefix>.gene_tpm.gct.gz` | Cohort-level gene TPM matrix in GCT format. |
| `count_gct` | `<prefix>.gene_reads.gct.gz` | Cohort-level gene read-count matrix in GCT format. |
| `exon_count_gct` | `<prefix>.exon_reads.gct.gz` | Exon-count matrix in GCT format when `merge_exons` is `true`. |
| `insert_size_hists` | `<prefix>.insert_size_hists.txt.gz` | Aggregated insert-size histograms when RNA-SeQC2 insert-size outputs are present. |
| `sample_count` | `sample_count.txt` | Number of samples validated from the manifest. |
| `batch_count` | `batch_count.txt` | Number of aggregation batches created. |

## Workflow and WDL

- Workflow: `rnaseqc2_aggregate_batched_workflow`
- Source WDL: [`rnaseqc2_aggregate_batched.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/GTEx/rnaseqc2_aggregate_batched.wdl)

## Processing stages

1. **Validate the manifest:** checks the manifest, determines the sample and batch counts, and records whether insert-size files are available.
2. **Aggregate batches:** downloads the RNA-SeQC2 outputs for each batch and merges TPM, read-count, metrics, and optional exon-count or insert-size files.
3. **Merge the cohort:** combines batch-level outputs into final cohort-level files using the requested prefix.

## Versioning

All `rnaseqc2_aggregate_batched_workflow` releases are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/GTEx/rnaseqc2_aggregate_batched.changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>