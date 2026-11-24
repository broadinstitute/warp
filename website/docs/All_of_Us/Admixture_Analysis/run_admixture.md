---
sidebar_position: 5
slug: /All_of_Us/Admixture_Analysis/run_admixture
title: Admixture Unsupervised
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.0.2](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/run_admixture.changelog.md) | January, 2026 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the run_admixture workflow

[`run_admixture`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/run_admixture.wdl) is a WDL workflow that runs [ADMIXTURE](http://dalexander.github.io/admixture/) in unsupervised mode to estimate global ancestry proportions from PLINK binary genotype data. It is the second step in the **Admixture Unsupervised** analysis path used in All of Us processing.

The workflow accepts PLINK binary files (`.bed/.bim/.fam`) produced by [`convert_vcf_to_plink_bed`](./convert_vcf_to_plink_bed.md) and runs ADMIXTURE to identify *K* ancestry components without reference population supervision. This approach allows for customization of the reference panel to better account for underrepresented populations compared to supervised methods.

ADMIXTURE outputs an ancestry proportion matrix (`.Q`) and a population allele frequency matrix (`.P`), which can be used for downstream population-level inference or construction of a pruned reference panel. For the ADMIXTURE manual and full parameter documentation, see the [ADMIXTURE documentation](http://dalexander.github.io/admixture/admixture-manual.pdf).

## Quickstart table

| Workflow Feature | Description | Source |
| :--: | :-- | :--: |
| Analysis type | Unsupervised admixture estimation | |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Input data type | PLINK binary (`.bed`, `.bim`, `.fam`) | |
| Output file format | `.Q` (ancestry proportions), `.P` (allele frequencies) | |
| Primary tool | ADMIXTURE | [ADMIXTURE](http://dalexander.github.io/admixture/) |
| Docker image | `mussmann/admixpipe:3.0` | |
| Part of analysis path | Admixture Unsupervised (Step 2 of 2) | |

## Set-up

### run_admixture installation and requirements

The workflow code can be downloaded by cloning the [WARP GitHub repository](https://github.com/broadinstitute/warp). For the latest release, please see the [run_admixture changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/run_admixture.changelog.md).

The pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant workflow management system.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `bed` | PLINK binary genotype file produced by [`convert_vcf_to_plink_bed`](./convert_vcf_to_plink_bed.md). | File |
| `bim` | PLINK variant information file. | File |
| `fam` | PLINK sample metadata file. | File |

The following values are currently set at the **task level** in the WDL and are not exposed as workflow inputs: `K_in` (default `6`), `num_cpus_in` (default `4`), and `mem_gb` (default `120`).

## run_admixture tasks and tools

This workflow calls a single task to perform unsupervised admixture clustering.

1. [Run ADMIXTURE unsupervised clustering](#1-run-admixture-unsupervised-clustering)

To see specific tool parameters, select the task WDL link in the table; then view the `command {}` section of the task in the WDL script.

| Task name and WDL link | Tool | Software | Description |
| --- | --- | --- | --- |
| [run_admixture](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/run_admixture.wdl) | ADMIXTURE | `/app/bin/admixture` (admixpipe:3.0) | Runs ADMIXTURE in unsupervised mode on PLINK binary files. Estimates *K* ancestry components using maximum likelihood. Outputs ancestry proportion and allele frequency matrices. |

### 1. Run ADMIXTURE unsupervised clustering

The task invokes ADMIXTURE on the input `.bed` file with the specified number of ancestry components (*K*) and CPU threads. ADMIXTURE uses maximum likelihood to estimate the proportion of each sample's ancestry attributable to each of the *K* inferred populations. Output files follow ADMIXTURE's standard naming convention based on the input file basename:

- `<basename>.<K>.Q` — ancestry proportion matrix (one row per sample, one column per ancestry component)
- `<basename>.<K>.P` — population allele frequency matrix (one row per variant, one column per ancestry component)

## Outputs

| Output variable name | Filename | Output format and description |
| --- | --- | --- |
| `admixture_Q` | `<basename>.<K>.Q` | Ancestry proportion matrix. Each row is a sample; each column is an inferred ancestry component. Values sum to 1 per row. |
| `admixture_P` | `<basename>.<K>.P` | Population allele frequency matrix. Each row is a variant; each column is the allele frequency in the corresponding inferred ancestry component. |

## Versioning

All `run_admixture` releases are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/run_admixture.changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
