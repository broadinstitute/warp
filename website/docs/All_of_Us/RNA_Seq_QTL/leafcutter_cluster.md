---
sidebar_position: 4
slug: /All_of_Us/RNA_Seq_QTL/leafcutter_cluster
title: Leafcutter Clustering
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.0.0](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/leafcutter_cluster.changelog.md) | July, 2025 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the Leafcutter Clustering workflow

[`leafcutter_cluster.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/leafcutter_cluster.wdl) defines workflow `leafcutter_cluster_workflow`, which clusters junction data and produces sQTL phenotype artifacts including bed/parquet matrices, phenotype groups, and PCs.

The workflow localizes junction files from a file-of-paths list, runs cluster preparation, and outputs leafcutter-ready downstream inputs.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `junc_files_list` | File containing one junction file path per line. | File |
| `exon_list` | Exon annotation list for clustering script. | File |
| `genes_gtf` | Gene GTF used for clustering and feature annotation. | File |
| `prefix` | Prefix used for all output filenames. | String |
| `sample_participant_lookup` | Sample-to-participant mapping file. | File |
| `min_clu_reads` | *(Optional)* minimum cluster reads threshold. | Int? |
| `min_clu_ratio` | *(Optional)* minimum cluster ratio threshold. | Float? |
| `max_intron_len` | *(Optional)* maximum intron length. | Int? |
| `num_pcs` | *(Optional)* number of PCs to compute. | Int? |
| `memory` | Runtime memory in GB. | Int |
| `disk_space` | Runtime disk size. | Int |
| `num_threads` | Runtime CPU threads. | Int |
| `num_preempt` | Runtime preemptible count. | Int |

## Outputs

| Output variable name | Filename, if applicable | Output format and description |
| --- | --- | --- |
| `counts_out` | `<prefix>_perind.counts.gz` | Per-individual intron cluster counts. |
| `counts_numers_out` | `<prefix>_perind_numers.counts.gz` | Per-individual numerator counts. |
| `clusters_pooled_out` | `<prefix>_pooled.gz` | Pooled cluster definitions. |
| `clusters_refined_out` | `<prefix>_refined.gz` | Refined cluster definitions. |
| `phenotype_groups_out` | `<prefix>.leafcutter.phenotype_groups.txt` | Phenotype group file for TensorQTL-style sQTL analysis. |
| `bed_parquet_out` | `<prefix>.leafcutter.bed.parquet` | Leafcutter phenotype matrix in parquet format. |
| `bed_out` | `<prefix>.leafcutter.bed.gz` | Leafcutter phenotype matrix in BED-like gzip format. |
| `bed_index_out` | `<prefix>.leafcutter.bed.gz.tbi` | Tabix index for `bed_out`. |
| `pcs_out` | `<prefix>.leafcutter.PCs.txt` | Phenotype PCs from clustering pipeline. |
| `leafcutter_pipeline_version` | `aou_9.0.0` | Workflow version string output. |

## Workflow and WDL

- Workflow: `leafcutter_cluster_workflow`
- Source WDL: [`leafcutter_cluster.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/leafcutter_cluster.wdl)

## Versioning

All `leafcutter_cluster` releases are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/leafcutter_cluster.changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
