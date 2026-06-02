---
sidebar_position: 6
slug: /All_of_Us/RNA_Seq_QTL/aggregate_susie_workflow
title: Aggregate SuSiE Workflow
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.0.2](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/AggregateSusieWorkflow.changelog.md) | January, 2026 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the Aggregate SuSiE workflow

[`AggregateSusieWorkflow.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/AggregateSusieWorkflow.wdl) merges per-phenotype SuSiE outputs and annotates the merged result with external reference resources.

It first aggregates all SuSiE parquet paths from a file-of-paths input, then annotates merged results using Gencode, allele-frequency, and functional annotation resources.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `SusieParquetsFOFN` | File-of-file-paths listing SuSiE parquet outputs to aggregate. | File |
| `Memory` | Runtime memory (GB). | Int |
| `OutputPrefix` | Prefix used for merged outputs. | String |
| `NumThreads` | CPU thread count for aggregation task. | Int |
| `GencodeGTF` | Gencode GTF for annotation context. | File |
| `PlinkAfreq` | Allele frequency table from PLINK processing. | File |
| `AnnotationPhyloP` | PhyloP annotation resource. | File |
| `AnnotationENCODE` | ENCODE cCRE annotation resource. | File |
| `AnnotationFANTOM5` | FANTOM5 annotation resource. | File |
| `AnnotationVEP` | VEP annotation table. | File |
| `AnnotationGnomad` | gnomAD annotation resource. | File |
| `AnnotationVEPIndex` | Index for VEP annotation table. | File |

## Outputs

| Output variable name | Filename, if applicable | Output format and description |
| --- | --- | --- |
| `AnnotatedMergedSusieParquet` | `<OutputPrefix>_SusieMerged.annotated.tsv` | Annotated merged SuSiE output table. |

## Workflow and WDL

- Workflow: `AggregateSusieWorkflow`
- Source WDL: [`AggregateSusieWorkflow.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/AggregateSusieWorkflow.wdl)

## Versioning

All `AggregateSusieWorkflow` releases are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/AggregateSusieWorkflow.changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
