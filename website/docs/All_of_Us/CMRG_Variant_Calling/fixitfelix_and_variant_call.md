---
sidebar_position: 2
slug: /All_of_Us/CMRG_Variant_Calling/fixitfelix_and_variant_call
title: FixItFelix and Variant Calling
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.0.1](https://github.com/broadinstitute/warp/blob/develop/all_of_us/cmrg/FixItFelixAndVariantCall.changelog.md) | January, 2026 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the FixItFelix and Variant Calling workflow

[`FixItFelixAndVariantCall`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/cmrg/FixItFelixAndVariantCall.wdl) performs CMRG-focused read extraction/remapping and variant calling on corrected alignments. It subsets the input alignment to true/false-duplication intervals, remaps reads against a masked GRCh38 reference with FixItFelix, and calls variants in true-location intervals.

The workflow supports either standard VCF output or optional gVCF output (via `generate_gvcf=true`).

## Quickstart table

| Pipeline Feature | Description | Source |
| :--: | :-- | :--: |
| Analysis type | CMRG remapping + variant calling | |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Genomic reference sequence | Masked GRCh38 (default inputs provided) | |
| Data input file format | CRAM/BAM + reference + interval BEDs | |
| Data output file format | VCF or gVCF (`.vcf.gz` / `.g.vcf.gz`) + index | |
| Primary software | GATK + FixItFelix + samtools + bwa | [GATK](https://gatk.broadinstitute.org/) |

## Set-up

### FixItFelix and Variant Calling installation and requirements

The workflow code can be downloaded by cloning the [WARP GitHub repository](https://github.com/broadinstitute/warp). For the latest release, please see the [FixItFelixAndVariantCall changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/cmrg/FixItFelixAndVariantCall.changelog.md).

The pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant workflow management system.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `cram_file` | Input alignment file (CRAM/BAM). | File |
| `cram_file_index` | Index for `cram_file` (`.crai`/`.bai`). | File |
| `original_ref_fasta` | Reference fasta used for the original alignment. | File |
| `original_ref_fasta_index` | FASTA index for original reference. | File |
| `original_ref_dict` | Sequence dictionary for original reference. | File |
| `masked_ref_fasta` | Masked GRCh38 reference for FixItFelix remapping. | File |
| `masked_ref_fasta_index` | FASTA index for masked reference. | File |
| `masked_ref_dict` | Dictionary for masked reference. | File |
| `masked_ref_amb` | BWA index `.amb` for masked reference. | File |
| `masked_ref_ann` | BWA index `.ann` for masked reference. | File |
| `masked_ref_bwt` | BWA index `.bwt` for masked reference. | File |
| `masked_ref_pac` | BWA index `.pac` for masked reference. | File |
| `masked_ref_sa` | BWA index `.sa` for masked reference. | File |
| `true_locations_intervals` | BED of true genomic locations used for final variant calling. | File |
| `false_duplications_intervals` | BED of false-duplication locations. | File |
| `combined_true_false_intervals` | BED combining true + false regions for read extraction/remapping. | File |
| `generate_gvcf` | Emit gVCF instead of VCF. Default: `false`. | Boolean |

## FixItFelix and Variant Calling tasks and tools

The workflow runs three tasks in sequence.

1. [Subset alignment to CMRG-relevant intervals](#1-subset-alignment-to-cmrg-relevant-intervals)
2. [Remap extracted reads with FixItFelix](#2-remap-extracted-reads-with-fixitfelix)
3. [Call variants on corrected BAM](#3-call-variants-on-corrected-bam)

| Task name and WDL link | Tool | Software | Description |
| --- | --- | --- | --- |
| [subset_cram](https://github.com/broadinstitute/warp/blob/develop/all_of_us/cmrg/FixItFelixAndVariantCall.wdl) | GATK PrintReads | `us.gcr.io/broad-gatk/gatk:4.4.0.0` | Subsets input alignment to combined true/false intervals and outputs BAM. |
| [FixItFelix](https://github.com/broadinstitute/warp/blob/develop/all_of_us/cmrg/FixItFelixAndVariantCall.wdl) | FixItFelix | `gcr.io/broad-dsde-methods/fixitfelix:1` | Extracts and remaps read pairs to masked GRCh38 reference. |
| [call_variants](https://github.com/broadinstitute/warp/blob/develop/all_of_us/cmrg/FixItFelixAndVariantCall.wdl) | GATK HaplotypeCaller + SelectVariants | `us.gcr.io/broad-gatk/gatk:4.4.0.0` | Calls variants in true-location intervals and optionally filters non-variant records for VCF mode. |

### 1. Subset alignment to CMRG-relevant intervals

`subset_cram` extracts reads overlapping `combined_true_false_intervals` to produce a smaller BAM for downstream remapping.

### 2. Remap extracted reads with FixItFelix

`FixItFelix` remaps extracted read pairs against the masked reference and outputs coordinate-sorted/indexed BAM.

### 3. Call variants on corrected BAM

`call_variants` runs HaplotypeCaller in DRAGEN mode over true-location intervals and emits either gVCF or filtered VCF.

## Outputs

| Output variable name | Filename, if applicable | Output format and description |
| --- | --- | --- |
| `output_vcf` | `filtered_<sample>.vcf.gz` or `<sample>.g.vcf.gz` | Final called variants (VCF or gVCF depending on `generate_gvcf`). |
| `output_vcf_index` | `filtered_<sample>.vcf.gz.tbi` or `<sample>.g.vcf.gz.tbi` | Tabix index for `output_vcf`. |
| `output_pipeline_version` | `aou_9.0.1` | Workflow version string output. |

## Versioning

All `FixItFelixAndVariantCall` releases are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/cmrg/FixItFelixAndVariantCall.changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
