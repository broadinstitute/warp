---
sidebar_position: 4
slug: /All_of_Us/Admixture_Analysis/convert_vcf_to_plink_bed
title: Admixture Unsupervised Preprocessing
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.0.1](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/convert_vcf_to_plink_bed.changelog.md) | November, 2025 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the convert_vcf_to_plink_bed workflow

[`convert_vcf_to_plink_bed`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/convert_vcf_to_plink_bed.wdl) is a WDL workflow that converts a merged VCF file into [PLINK binary format](https://www.cog-genomics.org/plink/1.9/formats#bed) (`.bed/.bim/.fam`). It is the first step in the **Admixture Unsupervised** analysis path used in All of Us processing.

The workflow accepts merged VCF shards produced by the ancestry pipeline and runs [PLINK](https://www.cog-genomics.org/plink/) to generate the binary genotype files required by the downstream [`run_admixture`](./run_admixture.md) workflow. PLINK is run with double-ID assignment and permissive chromosome handling to accommodate multi-ancestry cohort data.

This workflow must be run before [`run_admixture`](./run_admixture.md).

## Quickstart table

| Workflow Feature | Description | Source |
| :--: | :-- | :--: |
| Analysis type | VCF-to-PLINK format conversion | |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Input data type | Merged VCF + index | |
| Output file format | PLINK binary (`.bed`, `.bim`, `.fam`) | |
| Primary tool | PLINK | [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/) |
| Docker image | `mussmann/admixpipe:3.0` | |
| Part of analysis path | Admixture Unsupervised (Step 1 of 2) | |

## Set-up

### convert_vcf_to_plink_bed installation and requirements

The workflow code can be downloaded by cloning the [WARP GitHub repository](https://github.com/broadinstitute/warp). For the latest release, please see the [convert_vcf_to_plink_bed changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/convert_vcf_to_plink_bed.changelog.md).

The pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant workflow management system.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `prefix` | Output file prefix. Used as the base name for all PLINK output files. | String |
| `merged_vcf_shards` | Merged VCF file from the ancestry pipeline. | File |
| `merged_vcf_shards_idx` | Index file for the merged VCF. | File |

## convert_vcf_to_plink_bed tasks and tools

This workflow calls a single task to convert VCF input to PLINK binary format.

1. [Convert VCF to PLINK binary format](#1-convert-vcf-to-plink-binary-format)

To see specific tool parameters, select the task WDL link in the table; then view the `command {}` section of the task in the WDL script.

| Task name and WDL link | Tool | Software | Description |
| --- | --- | --- | --- |
| [convert_vcf_to_plink_bed](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/convert_vcf_to_plink_bed.wdl) | PLINK | `/app/bin/plink` (admixpipe:3.0) | Converts a merged VCF to PLINK binary format using `--make-bed`. Applies `--double-id` and `--allow-extra-chr` for compatibility with AoU cohort data. |

### 1. Convert VCF to PLINK binary format

The task invokes PLINK 1.9 with the following key flags:

| PLINK flag | Value | Notes |
| --- | --- | --- |
| `--vcf` | `<merged_vcf_shards>` | Input VCF file |
| `--make-bed` | — | Outputs binary PLINK format |
| `--double-id` | — | Sets both family ID and individual ID to the sample ID |
| `--allow-extra-chr` | — | Permits non-standard chromosome names |
| `--out` | `<prefix>` | Sets output file base name |

## Outputs

| Output variable name | Filename | Output format and description |
| --- | --- | --- |
| `bed` | `<prefix>.bed` | PLINK binary genotype file. |
| `bim` | `<prefix>.bim` | Variant information file (chromosome, variant ID, position, alleles). |
| `fam` | `<prefix>.fam` | Sample metadata file (family ID, individual ID, parental IDs, sex, phenotype). |

## Versioning

All `convert_vcf_to_plink_bed` releases are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/convert_vcf_to_plink_bed.changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
