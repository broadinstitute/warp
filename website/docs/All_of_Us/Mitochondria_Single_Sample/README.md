---
sidebar_position: 1
slug: /Pipelines/mitochondria_single_sample/README
---

# Mitochondria Single Sample Pipeline Overview

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.0.0](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/single_sample/mitochondria_single_sample.changelog.md) | May, 2026 | [WARP Pipelines](mailto:warp@broadinstitute.org) | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the Mitochondria Single Sample workflow

The Mitochondria Mitochondria Single Sample Pipeline ([`mitochondria_single_sample.wdl`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/single_sample/mitochondria_single_sample.wdl)) is a cloud-optimized WDL pipeline that processes per-sample mitochondrial DNA (mtDNA) data from hg38-aligned whole-genome sequencing (WGS) BAM or CRAM files. It produces high-quality per-sample mtDNA variant calls, base-level coverage metrics, and supporting QC statistics for use in downstream cohort-level analysis.

The pipeline is based on the [mtSwirl v2.5_MongoSwirl_Single](https://github.com/rahulg603/mtSwirl) methodology and implements a two-round alignment and variant calling strategy. In the first round, reads are aligned to the standard hg38 mitochondrial reference and variants are called with Mutect2. A personalized self-reference is then constructed from the called variants, and in the second round reads are realigned to this self-reference for improved variant calling accuracy. Variant calls from the second round are lifted back to hg38 coordinates for a consistent final output.

The pipeline optionally computes NUMT (nuclear mitochondrial DNA segment) coverage metrics for quality control purposes. Outputs from this pipeline — particularly the final VCF and per-base coverage metrics — serve as per-sample inputs to the [`MitochondriaMerge`](../MtCoverageMerge_Pipeline/README.md) cohort-level pipeline.

This pipeline was originally released as part of: Gupta, R., Kanai, M., Durham, T.J. et al. Nuclear genetic control of mtDNA copy number and heteroplasmy in humans. *Nature*, 2023. https://doi.org/10.1038/s41586-023-06426-5.

## Quickstart table

| Pipeline Feature | Description |
| :--: | :-- |
| Assay type | Whole-genome sequencing (mtDNA) |
| Overall workflow | Subsetting, alignment, variant calling, liftover, coverage metrics |
| Workflow language | WDL 1.0 |
| Genomic reference sequence | hg38 |
| Variant caller | Mutect2 (GATK) |
| Data input file format | BAM or CRAM (hg38-aligned) |
| Data output file format | VCF, TSV (coverage metrics) |

## Set-up

### Mitochondria Single Sample Pipeline installation and requirements

The pipeline code can be downloaded by cloning the [WARP GitHub repository](https://github.com/broadinstitute/warp). The WDL is located at [`all_of_us/mitochondria/mitochondria_single_sample.wdl`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/single_sample/mitochondria_single_sample.wdl).

The pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant workflow management system.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `wgs_aligned_input_bam_or_cram` | hg38-aligned WGS BAM or CRAM file | File |
| `wgs_aligned_input_bam_or_cram_index` | Index for the input BAM or CRAM | File? |
| `sample_name` | Sample identifier used in output filenames | String |
| `mt_interval_list` | Picard-style interval list for chrM (including putative NUMT intervals) | File |
| `nuc_interval_list` | Interval list for nuclear NUMT regions | File |
| `ref_fasta` | hg38 reference FASTA | File |
| `ref_fasta_index` | Index for the reference FASTA | File |
| `ref_dict` | Sequence dictionary for the reference FASTA | File |
| `mt_fasta` | Mitochondrial reference FASTA | File |
| `mt_fasta_index` | Index for the mitochondrial reference FASTA | File |
| `mt_dict` | Sequence dictionary for the mitochondrial reference FASTA | File |
| `blacklisted_sites` | BED file of artifact-prone sites to exclude | File |
| `blacklisted_sites_index` | Index for the blacklisted sites BED file | File |
| `control_region_shifted_reference_interval_list` | Interval list for the shifted control region reference | File |
| `non_control_region_interval_list` | Interval list for the non-control region | File |
| `HailLiftover` | Hail liftover script | File |
| `FaRenamingScript` | Script for FASTA renaming during self-reference construction | File |
| `CheckVariantBoundsScript` | Script to validate variant coordinates | File |
| `CheckHomOverlapScript` | Script to check homoplasmy overlaps | File |
| `JsonTools` | JSON utilities script | File |
| `force_manual_download` | Force manual download of input files | Boolean |
| `max_read_length` | Read length for optimization (default: 151) | Int? |
| `vaf_filter_threshold` | Hard threshold for filtering low-VAF sites | Float? |
| `f_score_beta` | F-score beta balancing recall vs. precision in filtering | Float? |
| `verifyBamID` | VerifyBamID contamination threshold | Float? |
| `compress_output_vcf` | Whether to compress output VCF | Boolean |
| `compute_numt_coverage` | Whether to compute NUMT coverage metrics | Boolean |
| `use_haplotype_caller_nucdna` | Whether to use HaplotypeCaller for nuclear DNA | Boolean |
| `skip_restore_hardclips` | Skip restoring hardclips (required for some CRAMs, e.g., AoU) | Boolean |
| `gatk_version` | GATK version string | String |
| `ucsc_docker` | Docker image for UCSC tools | String |
| `genomes_cloud_docker` | Docker image for genomes-in-the-cloud tools | String |
| `haplochecker_docker` | Docker image for haplochecker | String |
| `gatk_samtools_docker` | Docker image for GATK + samtools tasks | String |

## Mitochondria Single Sample Pipeline tasks and tools

The pipeline calls a series of sub-workflows and tasks. The following sections describe each major step:

1. [Subset and revert BAM to chrM](#1-subset-and-revert-bam-to-chrm)
2. [Round 1: align and call variants on standard reference](#2-round-1-align-and-call-variants-on-standard-reference)
3. [Construct self-reference](#3-construct-self-reference)
4. [Round 2: align and call variants on self-reference](#4-round-2-align-and-call-variants-on-self-reference)
5. [Liftover and collect outputs](#5-liftover-and-collect-outputs)
6. [(Optional) NUMT coverage](#6-optional-numt-coverage)

| Task / sub-workflow | Tool | Description |
| --- | --- | --- |
| [`MongoSubsetBamToChrMAndRevert`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/single_sample_subworkflows_and_tasks/MongoTasks_v2_5_Single.wdl) | GATK PrintReads, RevertSam | Subsets WGS BAM to chrM and NUMT intervals, reverts to unmapped BAM |
| [`AlignAndCallR1`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/single_sample_subworkflows_and_tasks/AlignAndCallR1_v2_5_Single.wdl) | BWA, Mutect2, HaplotypeCaller | Aligns to hg38 mtDNA reference, calls variants, estimates contamination and haplogroup |
| [`ProduceSelfReferenceFiles`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/single_sample_subworkflows_and_tasks/ProduceSelfReferenceFiles_v2_5_Single.wdl) | bcftools consensus, Picard, UCSC liftOver | Constructs a personalized self-reference FASTA and chain files from Round 1 variant calls |
| [`AlignAndCallR2`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/single_sample_subworkflows_and_tasks/AlignAndCallR2_v2_5_Single.wdl) | BWA, Mutect2 | Aligns to self-reference, calls variants, collects coverage metrics |
| [`MongoLiftoverVCFAndGetCoverage`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/single_sample_subworkflows_and_tasks/MongoTasks_v2_5_Single.wdl) | Picard LiftoverVcf | Lifts Round 2 VCF back to hg38 coordinates |
| [`MongoLiftoverSelfAndCollectOutputs`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/single_sample_subworkflows_and_tasks/MongoTasks_v2_5_Single.wdl) | Picard, custom scripts | Collects final coverage metrics in hg38 coordinates |
| [`NucCoverageAtEveryBase`](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/single_sample/mitochondria_single_sample.wdl) | Picard CollectHsMetrics | *(Optional)* Computes per-base NUMT coverage before and after realignment |

### 1. Subset and revert BAM to chrM

The `MongoSubsetBamToChrMAndRevert` task uses GATK PrintReads to subset the input WGS BAM or CRAM to the chrM and NUMT intervals, then reverts the subset to an unmapped BAM for realignment. Mean chrM coverage is computed at this step.

### 2. Round 1: align and call variants on standard reference

`AlignAndCallR1` aligns the unmapped chrM reads to the standard hg38 mitochondrial reference using BWA, marks duplicates, and calls SNPs and INDELs using Mutect2. HaplotypeCaller is optionally used for nuclear DNA variants. Contamination is estimated and a major haplogroup is assigned.

### 3. Construct self-reference

`ProduceSelfReferenceFiles` uses the Round 1 variant calls to construct a personalized self-reference FASTA via `bcftools consensus`. Chain files for lifting variant coordinates between the self-reference and hg38 are produced, along with updated interval lists.

### 4. Round 2: align and call variants on self-reference

`AlignAndCallR2` realigns the unmapped reads to the self-reference and calls variants again with Mutect2. Alignment to the self-reference reduces reference bias and improves genotype accuracy, particularly for heteroplasmic sites.

### 5. Liftover and collect outputs

`MongoLiftoverVCFAndGetCoverage` lifts the Round 2 VCF back to hg38 coordinates. `MongoLiftoverSelfAndCollectOutputs` collects per-base coverage metrics in hg38 space.

### 6. (Optional) NUMT coverage

When `compute_numt_coverage = true`, `NucCoverageAtEveryBase` uses Picard `CollectHsMetrics` to compute per-base NUMT coverage on the original reference, self-reference, and shifted self-reference for QC comparison.

## Outputs

| Output variable name | Description |
| --- | --- |
| `subset_bam` / `subset_bai` | chrM-subset BAM and index |
| `r1_vcf` / `r1_vcf_index` | Round 1 filtered VCF and index |
| `r1_nuc_vcf` / `r1_nuc_vcf_index` | Round 1 nuclear DNA VCF and index |
| `self_ref_vcf` / `self_ref_vcf_index` | Round 2 VCF called on self-reference |
| `final_vcf` | Final VCF lifted back to hg38 coordinates |
| `final_rejected_vcf` | Variants that failed liftover |
| `final_base_level_coverage_metrics` | Per-base mtDNA coverage in hg38 coordinates |
| `numt_base_level_coverage` | *(Optional)* Per-base NUMT coverage |
| `self_reference_fasta` | Personalized self-reference FASTA |
| `reference_to_self_ref_chain` | Chain file from hg38 → self-reference |
| `stats_outputs` | Summary statistics table |
| `mean_coverage` | Mean mtDNA coverage depth |
| `median_coverage` | Median mtDNA coverage depth |
| `major_haplogroup` | Assigned mtDNA haplogroup |
| `contamination` | Estimated contamination fraction |
| `mtdna_consensus_overlaps` | Number of overlapping consensus sites |

## Versioning and testing

All Mitochondria Single Sample Pipeline releases are documented in the [mitochondria_single_sample changelog](https://github.com/broadinstitute/warp/blob/master/all_of_us/mitochondria/single_sample/mitochondria_single_sample.changelog.md).

## Citing the Mitochondria Single Sample

When citing this pipeline, please cite the original publication:

Gupta, R., Kanai, M., Durham, T.J. et al. Nuclear genetic control of mtDNA copy number and heteroplasmy in humans. *Nature*, 2023. https://doi.org/10.1038/s41586-023-06426-5.

When citing WARP, please use:

Kylee Degatano, Aseel Awdeh, Robert Sidney Cox III, Wes Dingman, George Grant, Farzaneh Khajouei, Elizabeth Kiernan, Kishori Konwar, Kaylee L Mathews, Kevin Palis, Nikelle Petrillo, Geraldine Van der Auwera, Chengchen (Rex) Wang, Jessica Way. "Warp Analysis Research Pipelines: Cloud-optimized workflows for biological data processing and reproducible analysis." *Bioinformatics*, 2025; https://doi.org/10.1093/bioinformatics/btaf494.

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues).
