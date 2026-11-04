---
sidebar_position: 2
slug: /All_of_Us/RNA_Seq_QTL/rnaseq_aou
title: RNA-seq AoU Processing
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.0.1](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/GTEx/rnaseq_aou.changelog.md) | June, 2025 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the RNA-seq AoU Processing workflow

[`rnaseq_aou.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/GTEx/rnaseq_aou.wdl) defines workflow `rnaseq_pipeline_bam_workflow`, which orchestrates GTEx-style RNA-seq processing from aligned reads through expression quantification and QC task chains.

It calls `samtofastq`, `star`, `rsem`, `markduplicates`, and `rnaseqc2` component WDLs to produce alignment/quantification artifacts used in downstream eQTL/sQTL workflows.

## Quickstart table

| Pipeline Feature | Description | Source |
| :--: | :-- | :--: |
| Analysis type | RNA alignment and quantification orchestration | |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Data input file format | task-level inputs provided through imported component workflows | |
| Data output file format | task-level outputs from STAR/RSEM/RNA-SeQC components | |
| Primary software | STAR, RSEM, RNA-SeQC (via imported WDLs) | |

## Set-up

### RNA-seq AoU Processing installation and requirements

The workflow code can be downloaded by cloning the [WARP GitHub repository](https://github.com/broadinstitute/warp). For the latest release, please see the [rnaseq_aou changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/GTEx/rnaseq_aou.changelog.md).

The pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant workflow management system.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `prefix` | Prefix passed to component RNA-seq tasks for naming and sample identity propagation. | String |

## RNA-seq AoU Processing tasks and tools

This workflow is an orchestration wrapper around imported GTEx-derived WDLs.

1. [Convert alignment to FASTQ](#1-convert-alignment-to-fastq)
2. [Align and quantify transcripts](#2-align-and-quantify-transcripts)
3. [Mark duplicates and run RNA QC](#3-mark-duplicates-and-run-rna-qc)

| Task name and WDL link | Tool | Software | Description |
| --- | --- | --- | --- |
| [samtofastq](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/GTEx/rnaseq_aou.wdl) | imported task | see `samtofastq.wdl` | Generates FASTQ from input alignment. |
| [star](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/GTEx/rnaseq_aou.wdl) | imported task | see `star.wdl` | Performs STAR alignment. |
| [rsem](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/GTEx/rnaseq_aou.wdl) | imported task | see `rsem.wdl` | Generates expression quantification outputs. |
| [markduplicates](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/GTEx/rnaseq_aou.wdl) | imported task | see `markduplicates.wdl` | Marks duplicate reads in aligned BAM. |
| [rnaseqc2](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/GTEx/rnaseq_aou.wdl) | imported task | see `rnaseqc2.wdl` | Produces RNA quality metrics. |

### 1. Convert alignment to FASTQ

`samtofastq` emits paired FASTQs for downstream alignment.

### 2. Align and quantify transcripts

`star` and `rsem` generate alignments and expression quantifications.

### 3. Mark duplicates and run RNA QC

`markduplicates` and `rnaseqc2` produce cleaned BAMs and QC outputs used by downstream analyses.

## Outputs

This wrapper workflow does not declare a top-level `output` block in the WDL. Outputs are available from called component tasks in execution metadata.

## Versioning

All `rnaseq_aou` releases are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/GTEx/rnaseq_aou.changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
