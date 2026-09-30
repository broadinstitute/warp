---
sidebar_position: 2
slug: /All_of_Us/Ancestry_Analysis/determine_hq_sites
title: Determine HQ Sites
className: aou-doc-page
---

<div className="aou-folder-text">

## Introduction to the Determine HQ Sites workflow

[`determine_hq_sites`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/determine_hq_sites.wdl) is a WDL workflow for creating a reusable set of high-quality variant sites for ancestry analysis. It is intended to be run infrequently on training VCFs, such as HGDP and 1000 Genomes data, before the downstream ancestry workflows are run.

The workflow processes VCF shards in parallel, filters them to passing biallelic SNPs with allele frequency greater than 0.1% and call rate greater than 99%, applies the requested intervals, performs LD pruning, and merges the resulting shards. It produces both sites-only and full-genotype VCFs, along with interval files for downstream filtering.

## Quickstart table

| Pipeline Feature | Description | Source |
| :--- | :--- | :--- |
| Analysis type | High-quality training-site selection and LD pruning | |
| Workflow language | WDL 1.0 | [openWDL](https://openwdl.org/) |
| Genomic reference sequence | GRCh38 / hg38 | |
| Data input file format | BGZF-compressed VCF shards and a drop-annotation arguments file | |
| Data output file format | BGZF-compressed VCFs, VCF indexes, interval list, and UCSC BED | |
| Primary software | GATK, Hail, bcftools | [GATK](https://gatk.broadinstitute.org/), [Hail](https://hail.is/), [bcftools](https://samtools.github.io/bcftools/) |

## Set-up

The workflow code can be downloaded by cloning the [WARP GitHub repository](https://github.com/broadinstitute/warp). It can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant workflow management system.

This workflow assumes that the input VCFs are chromosome-sharded and supplied in sorted order. The input data should use the GRCh38/hg38 reference. Because LD pruning runs on a large single machine for each shard, this workflow may require substantial compute resources.

## Inputs

| Input variable name | Description | Type |
| --- | --- | --- |
| `vcfs` | Ordered training VCF shards to filter and merge. | Array[File] |
| `drop_info_annotations_param_file` | GATK arguments file listing INFO annotations to remove from the output VCFs. | File |
| `output_prefix` | Prefix used for the merged sites-only output and related files. | String |
| `intervals` | Intervals to consider when selecting high-quality sites, such as calling intervals. | File |
| `intersecting_intervals` | *(Optional)* Additional intervals, such as exome targets, used to further restrict the selected sites. | File? |

## Workflow steps

1. **Filter each VCF shard.** GATK `SelectVariants` retains passing biallelic SNPs with `AF > 0.001` and no more than 1% missing genotypes, optionally restricted to an additional interval file.
2. **LD-prune each filtered shard.** Hail imports each VCF using the GRCh38 reference and retains variants selected by `hl.ld_prune` with `r2 = 0.1` and a 500,000 bp window. Shards with zero or one variant bypass pruning.
3. **Create indexes and merge outputs.** The workflow indexes the sites-only and full-genotype outputs, merges the sites-only VCFs with GATK, concatenates the full-genotype BGZF VCFs with bcftools, and creates an interval list and UCSC BED file from the merged sites-only VCF.

## Outputs

| Output variable name | Filename, if applicable | Output format and description |
| --- | --- | --- |
| `vcf_bgz_merged_sites_only` | `<output_prefix>.vcf.gz` | Merged sites-only VCF containing the LD-pruned high-quality sites. |
| `vcf_bgz_merged_sites_only_index` | `<output_prefix>.vcf.gz.tbi` | Index for the merged sites-only VCF. |
| `vcf_bgz_merged_full` | `<output_prefix>_full.vcf.bgz` | Merged full-genotype VCF containing the LD-pruned sites. |
| `vcf_bgz_merged_full_index` | `<output_prefix>_full.vcf.bgz.tbi` | Index for the merged full-genotype VCF. |
| `vcf_bgz_fulls` | Per-input-shard `.full.vcf.bgz` files | Array of per-shard full-genotype VCFs before the final concatenation. |
| `vcf_bgz_fulls_index` | Per-input-shard `.full.vcf.bgz.tbi` files | Array of indexes for the per-shard full-genotype VCFs. |
| `hq_site_ucsc_bed` | `<output_prefix>.vcf.gz.ucsc.bed` | UCSC BED representation of the merged high-quality sites. |
| `hq_site_interval_list` | `<output_prefix>.vcf.gz.interval_list` | Picard/GATK interval-list representation of the merged high-quality sites. |

## Relationship to downstream workflows

This workflow creates the training high-quality site set. Use its sites-only and full-genotype outputs as inputs to the ancestry preparation workflow, such as [`determine_hq_sites_intersection`](./determine_hq_sites_intersection.md), which intersects the training sites with the target dataset before [`run_ancestry`](./run_ancestry.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>