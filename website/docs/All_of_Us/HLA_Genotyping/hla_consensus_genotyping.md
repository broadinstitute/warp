---
sidebar_position: 2
slug: /All_of_Us/HLA_Genotyping/hla_consensus_genotyping
title: HLA Consensus Genotyping
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.0.0](https://github.com/broadinstitute/warp/blob/develop/all_of_us/hla/HLAGenotyping_changelog.md) | August, 2025 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the HLA Consensus Genotyping workflow

[`HLAGenotyping`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/hla/HLAGenotyping.wdl) performs HLA typing from aligned sequencing reads. It extracts HLA-region reads, runs HLA-HD, and conditionally runs Polysolver/OptiType when two-field uncertainty is detected, then builds a consensus callset.

The workflow is designed for GRCh38-aligned data and outputs harmonized consensus calls suitable for downstream aggregation.

## Quickstart table

| Pipeline Feature | Description | Source |
| :--: | :-- | :--: |
| Analysis type | Per-sample HLA genotyping and consensus calling | |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Data input file format | BAM/CRAM + reference + HLA intervals + helper scripts | |
| Data output file format | TSV-style HLA result tables | |
| Primary software | GATK, HLA-HD, Polysolver, OptiType, samtools | [GATK](https://gatk.broadinstitute.org/) |

## Set-up

### HLA Consensus Genotyping installation and requirements

The workflow code can be downloaded by cloning the [WARP GitHub repository](https://github.com/broadinstitute/warp). For the latest release, please see the [HLAGenotyping changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/hla/HLAGenotyping_changelog.md).

The pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant workflow management system.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `gatk_docker` | Docker image containing GATK/samtools tools. | String |
| `hlahd_docker` | Docker image containing HLA-HD. | String |
| `polysolver_docker` | Docker image containing Polysolver. | String |
| `optitype_docker` | Docker image containing OptiType. | String |
| `original_bam` | Input BAM or CRAM. | File |
| `original_bam_idx` | Index for `original_bam`. | File |
| `ref_fasta` | Reference FASTA. | File |
| `ref_fai` | FASTA index. | File |
| `ref_dict` | Reference dictionary. | File |
| `hla_intervals` | Interval list for HLA region extraction. | File |
| `convert_alleles_python_script` | Helper script for allele conversion/normalization. | File |
| `count_two_field_alleles_python_script` | Helper script to count two-field calls in HLA-HD output. | File |
| `hla_groups_file` | HLA group mapping file used during conversion. | File |
| `gcs_project_for_requester_pays` | *(Optional)* requester-pays project for cloud reads. | String? |
| `EMPTY_STRING_HACK` | *(Optional)* Terra compatibility helper for optional empty values. | File? |

## HLA Consensus Genotyping tasks and tools

The workflow chains extraction, primary typing, conditional fallback typing, and consensus.

1. [Extract HLA-only reads and FASTQs](#1-extract-hla-only-reads-and-fastqs)
2. [Run HLA-HD typing](#2-run-hla-hd-typing)
3. [Run conditional Polysolver/OptiType and consensus](#3-run-conditional-polysolveroptitype-and-consensus)

| Task name and WDL link | Tool | Software | Description |
| --- | --- | --- | --- |
| [MakeHLAOnlyBamsAndFastqs](https://github.com/broadinstitute/warp/blob/develop/all_of_us/hla/HLAGenotyping.wdl) | GATK + samtools | user-supplied `gatk_docker` | Extracts HLA intervals, sorts/indexes HLA BAM, and emits paired FASTQs. |
| [HLAHD](https://github.com/broadinstitute/warp/blob/develop/all_of_us/hla/HLAGenotyping.wdl) | HLA-HD + Python post-processing | user-supplied `hlahd_docker` | Performs primary HLA typing and conversion to harmonized output format. |
| [Polysolver](https://github.com/broadinstitute/warp/blob/develop/all_of_us/hla/HLAGenotyping.wdl) | Polysolver | user-supplied `polysolver_docker` | Conditional secondary caller for low-confidence loci. |
| [Optitype](https://github.com/broadinstitute/warp/blob/develop/all_of_us/hla/HLAGenotyping.wdl) | OptiType | user-supplied `optitype_docker` | Conditional secondary caller for A/B/C loci. |
| [Consensus](https://github.com/broadinstitute/warp/blob/develop/all_of_us/hla/HLAGenotyping.wdl) | custom merge logic | workflow-internal | Produces final consensus callset across callers. |

### 1. Extract HLA-only reads and FASTQs

`MakeHLAOnlyBamsAndFastqs` slices input alignment to HLA intervals and prepares paired FASTQ inputs for downstream callers.

### 2. Run HLA-HD typing

`HLAHD` runs primary typing and outputs converted/harmonized calls plus two-field count for conditional branching.

### 3. Run conditional Polysolver/OptiType and consensus

When two-field calls are present, Polysolver and OptiType are run; `Consensus` combines outputs with HLA-HD to produce final consensus results.

## Outputs

| Output variable name | Filename, if applicable | Output format and description |
| --- | --- | --- |
| `hlahd_raw_result` | task-produced text file | Raw normalized HLA-HD result table. |
| `hlahd_converted_result` | task-produced text file | Converted HLA-HD calls in harmonized allele naming. |
| `hlahd_two_field_count` | integer value | Count of two-field loci from HLA-HD output. |
| `hlahd_overriden` | integer value | Number of loci overridden during consensus generation. |
| `consensus` | consensus text file | Final per-sample consensus HLA calls. |
| `optitype_result` | optional text file | OptiType output when conditional branch runs. |
| `polysolver_result` | optional text file | Polysolver output when conditional branch runs. |

## Versioning

All `HLAGenotyping` releases are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/hla/HLAGenotyping_changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
