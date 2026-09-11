---
sidebar_position: 3
slug: /All_of_Us/RNA_Seq_QTL/leafcutter_bam_to_junc
title: Leafcutter BAM to Junctions
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.0.0](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/leafcutter_bam_to_junc.changelog.md) | July, 2025 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the Leafcutter BAM to Junctions workflow

[`leafcutter_bam_to_junc.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/leafcutter_bam_to_junc.wdl) defines workflow `leafcutter_bam_to_junc_workflow`, which extracts splice junctions from BAM inputs for sQTL processing.

Reads are filtered for mapping quality and WASP tags, then `regtools junctions extract` produces compressed junction outputs per sample.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `bam_file` | Input aligned BAM file. | File |
| `sample_id` | Sample identifier used in output naming. | String |
| `strand_specificity` | *(Optional)* strand-specificity setting for regtools (`0` default). | Int? |
| `user_defined_threads` | *(Optional)* thread override (default `1`). | Int? |
| `user_defined_preempt` | *(Optional)* preemptible count override (default `1`). | Int? |

## Outputs

| Output variable name | Filename, if applicable | Output format and description |
| --- | --- | --- |
| `junc_file_final` | `<sample_id>.regtools_junc.txt.gz` | Gzipped splice junction file extracted from filtered BAM reads. |
| `pipeline_version_out` | `aou_9.0.0` | Workflow version string output. |

## Workflow and WDL

- Workflow: `leafcutter_bam_to_junc_workflow`
- Source WDL: [`leafcutter_bam_to_junc.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/leafcutter_bam_to_junc.wdl)

## Versioning

All `leafcutter_bam_to_junc` releases are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/leafcutter_bam_to_junc.changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
