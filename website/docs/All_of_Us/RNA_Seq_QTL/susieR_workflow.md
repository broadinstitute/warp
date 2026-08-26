---
sidebar_position: 5
slug: /All_of_Us/RNA_Seq_QTL/susieR_workflow
title: SuSiE Fine-Mapping Workflow
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.1.0](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/susieR_workflow.changelog.md) | August, 2026 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the SuSiE Fine-Mapping workflow

[`susieR_workflow.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/susieR_workflow.wdl) performs cis-window fine-mapping using SuSiE from TensorQTL outputs and dosage matrices.

Given an `InputMappingTsv` (one row per ancestry) and a target `Ancestry`, the workflow first looks up the ancestry-specific genotype dosage, covariate, sample list, phenotype BED, and TensorQTL permutation file paths, subsets those phenotype/genotype inputs for a target phenotype ID, then runs an R-based SuSiE script to produce parquet outputs including credible set summaries and full fine-mapping tables.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `InputMappingTsv` | Tab-separated file, one row per ancestry, providing `gs://` paths for `GenotypeDosage`, `GenotypeDosagei`, `PhenotypePCsOut`, `PlinkAF`, `QtlCovariates`, `Sample_list`, `VCF`, `cis_qtl`, `genotype_pcs`, `pgen`, `phenotype_bed`, and `psam`, keyed by a `mapping_inputs_id` ancestry code (AFR, AMR, COMB, EUR, SAS, EAS, MID). | File |
| `Ancestry` | Ancestry code used to select the row of `InputMappingTsv` to use for this run. | String |
| `CisDistance` | Cis-window distance parameter. | Int |
| `susie_rscript` | SuSiE runner script path. | File |
| `memory` | Runtime memory (GB). | Int |
| `NumPrempt` | Runtime preemptible count. | Int |
| `OutputPrefix` | General output prefix. | String |
| `PhenotypeID` | Target phenotype ID for subsetting and fine-mapping. | String |

## Outputs

| Output variable name | Filename, if applicable | Output format and description |
| --- | --- | --- |
| `SusieParquet` | `<PhenotypeID>.parquet` | Primary SuSiE fine-mapping parquet output. |
| `SusielbfParquet` | `<PhenotypeID>.lbf_variable.parquet` | Variant-level log-Bayes-factor parquet output. |
| `FullSusieParquet` | `<PhenotypeID>.full_susie.parquet` | Expanded full SuSiE results parquet. |
| `SubsetBed` | `<PhenotypeID>.bed.bgz` | Phenotype-specific subset bed used during run. |
| `SubsetDosages` | `<PhenotypeID>.dose.tsv.gz` | Dosage matrix subset for selected phenotype interval. |
| `SubsetDosagesIndex` | `<PhenotypeID>.dose.tsv.gz.tbi` | Tabix index for `SubsetDosages`. |

## Workflow and WDL

- Workflow: `susieR_workflow`
- Source WDL: [`susieR_workflow.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/susieR_workflow.wdl)

## Versioning

All `susieR_workflow` releases are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/susieR_workflow.changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
