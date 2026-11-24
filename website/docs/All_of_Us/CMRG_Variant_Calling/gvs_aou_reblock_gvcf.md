---
sidebar_position: 3
slug: /All_of_Us/CMRG_Variant_Calling/gvs_aou_reblock_gvcf
title: GVS AoU Reblock gVCF
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.0.0](https://github.com/broadinstitute/warp/blob/develop/all_of_us/cmrg/GvsAoUReblockGvcf.changelog.md) | August, 2025 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the GVS AoU Reblock gVCF workflow

[`GvsAoUReblockGvcf`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/cmrg/GvsAoUReblockGvcf.wdl) reblocks existing per-sample gVCFs using GATK `ReblockGVCF` and optionally copies results to site-specific AoU research buckets.

This workflow is commonly used as a preparation step before external joint calling systems (e.g., GVS-based aggregation workflows).

## Quickstart table

| Pipeline Feature | Description | Source |
| :--: | :-- | :--: |
| Analysis type | gVCF reblocking and optional site-bucket copy | |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Data input file format | gVCF path + reference bundle | |
| Data output file format | Reblocked gVCF + index | |
| Primary software | GATK ReblockGVCF + gsutil | [GATK](https://gatk.broadinstitute.org/) |

## Set-up

### GVS AoU Reblock gVCF installation and requirements

The workflow code can be downloaded by cloning the [WARP GitHub repository](https://github.com/broadinstitute/warp). For the latest release, please see the [GvsAoUReblockGvcf changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/cmrg/GvsAoUReblockGvcf.changelog.md).

The pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant workflow management system.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `gvcf` | GCS path to input gVCF (`.g.vcf.gz`). | String |
| `gvcf_index` | *(Optional)* Index path for input gVCF. Defaults to `gvcf + ".tbi"`. | String? |
| `ref_dict` | Reference sequence dictionary. | File |
| `ref_fasta` | Reference FASTA. | File |
| `ref_fasta_index` | FASTA index. | File |
| `requester_pays_project` | *(Optional)* GCS requester-pays project passed to GATK. | String? |
| `site_id` | *(Optional)* Destination site code (`bi`, `bcm`, `uw`). | String? |
| `docker_image` | Docker image for GATK ReblockGVCF. Default: `us.gcr.io/broad-gatk/gatk:4.2.6.1`. | String |

## GVS AoU Reblock gVCF tasks and tools

This workflow runs a single task.

1. [Reblock and optionally copy output](#1-reblock-and-optionally-copy-output)

| Task name and WDL link | Tool | Software | Description |
| --- | --- | --- | --- |
| [ReblockAndCopy](https://github.com/broadinstitute/warp/blob/develop/all_of_us/cmrg/GvsAoUReblockGvcf.wdl) | GATK ReblockGVCF | configurable (`docker_image`) | Reblocks input gVCF, then optionally copies output and index to a site-specific bucket prefix. |

### 1. Reblock and optionally copy output

`ReblockAndCopy` runs GATK with fixed GQ bands and optional requester-pays settings. If `site_id` is provided, output files are copied to the mapped destination bucket path.

## Outputs

| Output variable name | Filename, if applicable | Output format and description |
| --- | --- | --- |
| `reblocked_gvcf` | `<basename>.reblocked.g.vcf.gz` | Reblocked gVCF output (optionally delocalized to site-specific bucket path). |
| `reblocked_gvcf_index` | `<basename>.reblocked.g.vcf.gz.tbi` | Tabix index for reblocked gVCF. |

## Versioning

All `GvsAoUReblockGvcf` releases are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/cmrg/GvsAoUReblockGvcf.changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
