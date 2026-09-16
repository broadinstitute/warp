---
sidebar_position: 2
slug: /All_of_Us/Ancestry_Analysis/vds_to_vcf
title: VDS to VCF
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.0.0](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/vds_to_vcf.changelog.md) | July, 2025 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the VDS to VCF workflow

[`vds_to_vcf`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/vds_to_vcf.wdl) is a WDL workflow that converts a Hail Variant Dataset (VDS) into per-contig VCF outputs for downstream ancestry analysis. It is designed for large All of Us callsets aligned to GRCh38 and supports scatter-based chromosome processing for scalability.

The workflow repartitions the VDS, filters by contig and BED intervals, densifies to a matrix table, and exports both full and sites-only VCFs with Tabix indexes. It also writes two file-of-filenames (FOFN) manifests listing the generated VCFs and index files for downstream workflows.

## Quickstart table

| Pipeline Feature | Description | Source |
| :--: | :-- | :--: |
| Analysis type | VDS conversion and interval filtering for ancestry preprocessing | |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Genomic reference sequence | GRCh38 | |
| Data input file format | VDS + BED + contig list | |
| Data output file format | VCF BGZF + Tabix index + FOFN text files | |
| Primary software | Hail, GATK, bcftools | [Hail](https://hail.is/), [GATK](https://gatk.broadinstitute.org/), [bcftools](https://samtools.github.io/bcftools/) |

## Set-up

### VDS to VCF installation and requirements

The workflow code can be downloaded by cloning the [WARP GitHub repository](https://github.com/broadinstitute/warp). For the latest release, please see the [vds_to_vcf changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/vds_to_vcf.changelog.md).

The pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant workflow management system.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `vds_gs_url` | Google Cloud Storage path to the input Hail Variant Dataset. | String |
| `bed_gs_url` | Google Cloud Storage path to BED intervals used for filtering. | String |
| `n_partitions` | Number of partitions to apply to the VDS before processing. Default: `2000`. | Int |
| `output_prefix` | Prefix used for all generated output files. | String |
| `contigs` | Ordered list of contigs/chromosomes to process in scatter mode. | Array[String] |

## VDS to VCF tasks and tools

The VDS to VCF workflow calls two tasks to process each contig and create output manifests.

1. [Process VDS per contig](#1-process-vds-per-contig)
2. [Create output file manifests](#2-create-output-file-manifests)

To see specific tool parameters, select the task WDL link in the table; then view the `command {}` section of the task in the WDL script.

| Task name and WDL link | Tool | Software | Description |
| --- | --- | --- | --- |
| [process_vds](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/vds_to_vcf.wdl) | Hail, Python | `hailgenetics/hail:0.2.134-py3.11` | Repartitions VDS, filters by chromosome and BED intervals, and exports full + sites-only VCFs with indexes. |
| [create_fofn](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/vds_to_vcf.wdl) | Shell | `us.gcr.io/broad-gatk/gatk:4.2.6.1` | Writes text manifests listing output VCF and index file paths. |

### 1. Process VDS per contig

For each contig in `contigs`, the workflow calls `process_vds` to filter and export outputs named with `<output_prefix>.<contig>`. Each shard generates a full VCF and a sites-only VCF, each with a `.tbi` index.

### 2. Create output file manifests

After scatter completion, `create_fofn` writes two text files (`.fofn1.txt` and `.fofn2.txt`) containing the list of full VCFs and VCF index files.

## Outputs

| Output variable name | Filename, if applicable | Output format and description |
| --- | --- | --- |
| `vcfs` | `<output_prefix>.<contig>.vcf.bgz` | Array of per-contig full VCF files. |
| `vcfs_tbis` | `<output_prefix>.<contig>.vcf.bgz.tbi` | Array of Tabix indexes for full VCF files. |
| `vcfs_so` | `<output_prefix>.<contig>.so.vcf.bgz` | Array of per-contig sites-only VCF files. |
| `vcfs_so_tbis` | `<output_prefix>.<contig>.so.vcf.bgz.tbi` | Array of Tabix indexes for sites-only VCF files. |
| `vcfs_list` | `<output_prefix>.fofn1.txt` | Text file listing full VCF output paths. |
| `vcfs_idx_list` | `<output_prefix>.fofn2.txt` | Text file listing full VCF index output paths. |

## Versioning

All `vds_to_vcf` releases are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/vds_to_vcf.changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
