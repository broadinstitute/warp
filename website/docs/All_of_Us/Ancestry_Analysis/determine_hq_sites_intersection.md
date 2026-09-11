---
sidebar_position: 3
slug: /All_of_Us/Ancestry_Analysis/determine_hq_sites_intersection
title: Determine HQ Sites Intersection
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.0.0](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/determine_hq_sites_intersection.changelog.md) | September, 2025 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the Determine HQ Sites Intersection workflow

[`determine_hq_sites_intersection`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/determine_hq_sites_intersection.wdl) is a WDL workflow that identifies high-quality variant sites shared between a training callset and a target dataset, then filters both datasets to the shared intersection. This ensures downstream ancestry inference runs on a consistent set of informative sites.

The workflow processes input VCF shards, applies SNP and quality filters, intersects each shard with training sites-only variants, merges intersections, and generates three final outputs: the shared HQ sites-only VCF, filtered merged input-data VCF, and filtered training-set VCF.

## Quickstart table

| Pipeline Feature | Description | Source |
| :--: | :-- | :--: |
| Analysis type | Shared high-quality site selection and filtering | |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Data input file format | VCF BGZF + Tabix index (training + data shards) | |
| Data output file format | Filtered VCF BGZF + Tabix index | |
| Primary software | GATK, bcftools | [GATK](https://gatk.broadinstitute.org/), [bcftools](https://samtools.github.io/bcftools/) |

## Set-up

### Determine HQ Sites Intersection installation and requirements

The workflow code can be downloaded by cloning the [WARP GitHub repository](https://github.com/broadinstitute/warp). For the latest release, please see the [determine_hq_sites_intersection changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/determine_hq_sites_intersection.changelog.md).

The pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant workflow management system.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `training_vcf_bgz` | Full training VCF (with sample genotypes). | File |
| `training_vcf_bgz_idx` | Index for `training_vcf_bgz`. | File |
| `training_vcf_so_bgz` | Sites-only training VCF corresponding to `training_vcf_bgz`. | File |
| `training_vcf_so_bgz_idx` | Index for `training_vcf_so_bgz`. | File |
| `ordered_vcf_shards_in` | Ordered list of input dataset VCF shards. | Array[File] |
| `ordered_vcf_shards_idx_in` | Ordered list of VCF indexes corresponding to `ordered_vcf_shards_in`. | Array[File] |
| `ordered_vcf_shards_list` | *(Optional)* Text file listing shard VCF paths; overrides `ordered_vcf_shards_in` when provided. | File? |
| `ordered_vcf_shards_idx_list` | *(Optional)* Text file listing shard index paths; overrides `ordered_vcf_shards_idx_in` when provided. | File? |
| `final_output_prefix` | Output prefix for workflow-level naming. | String |
| `service_account_json` | *(Optional)* Service account key path used for requester-protected data localization. | String? |
| `intersecting_intervals` | *(Optional)* Additional intervals (e.g., exome targets) for intersection filtering. | File? |

## Determine HQ Sites Intersection tasks and tools

The workflow performs shard-level filtering/intersection followed by merge/filter steps to produce shared HQ callsets.

1. [Filter and sites-only transform each data shard](#1-filter-and-sites-only-transform-each-data-shard)
2. [Intersect shards with training HQ sites](#2-intersect-shards-with-training-hq-sites)
3. [Merge intersection sites and filter datasets](#3-merge-intersection-sites-and-filter-datasets)

To see specific tool parameters, select the task WDL link in the table; then view the `command {}` section of the task in the WDL script.

| Task name and WDL link | Tool | Software | Description |
| --- | --- | --- | --- |
| [sitesOnlyAndHQFilterVcf](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/determine_hq_sites_intersection.wdl) | GATK SelectVariants | `us.gcr.io/broad-gatk/gatk:4.2.0.0` | Filters each input shard to high-quality biallelic SNPs and writes a sites-only VCF. |
| [intersect_vcfs_as_sites_only](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/determine_hq_sites_intersection.wdl) | GATK SelectVariants | `us.gcr.io/broad-gatk/gatk:4.2.0.0` | Intersects each filtered shard with the training sites-only VCF. |
| [merge_vcf_bgzs](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/determine_hq_sites_intersection.wdl) | bcftools | `mgibio/bcftools-cwl:1.12` | Concatenates and sorts VCF shards into a merged VCF. |
| [filter_by_sites_only](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/determine_hq_sites_intersection.wdl) | GATK SelectVariants | `us.gcr.io/broad-gatk/gatk:4.2.0.0` | Filters full VCFs to the final intersected sites-only set. |

### 1. Filter and sites-only transform each data shard

Each shard is filtered for biallelic SNPs, allele frequency, missingness, and optional interval overlap, then exported as a sites-only VCF.

### 2. Intersect shards with training HQ sites

Each filtered shard is intersected against `training_vcf_so_bgz` to keep only sites shared with the training HQ set.

### 3. Merge intersection sites and filter datasets

Intersected shard sites are merged into a unified sites-only VCF. That merged site set is then used to filter both the training full VCF and all input data shards; filtered shards are merged into one data VCF.

## Outputs

| Output variable name | Filename, if applicable | Output format and description |
| --- | --- | --- |
| `hq_variants_intersection` | `merged_sites_only_intersection.vcf.bgz` | Sites-only VCF of shared high-quality variants between training and input datasets. |
| `hq_variants_intersection_idx` | `merged_sites_only_intersection.vcf.bgz.tbi` | Tabix index for `hq_variants_intersection`. |
| `merged_vcf_shards` | `merged_data_shards.vcf.bgz` | Merged input-data VCF filtered to shared HQ sites. |
| `merged_vcf_shards_idx` | `merged_data_shards.vcf.bgz.tbi` | Tabix index for `merged_vcf_shards`. |
| `filtered_training_set` | `full_training_sites_filtered.0.vcf.bgz` | Training full VCF filtered to shared HQ sites. |
| `filtered_training_set_idx` | `full_training_sites_filtered.0.vcf.bgz.tbi` | Tabix index for `filtered_training_set`. |

## Versioning

All `determine_hq_sites_intersection` releases are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/determine_hq_sites_intersection.changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
