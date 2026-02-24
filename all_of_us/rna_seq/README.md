# All of Us RNA-seq eQTL and sQTL Analysis Pipeline

This README describes the end-to-end workflow for preparing genotypes, generating RNA expression and splicing phenotypes, computing covariates, running cis-QTL analysis with TensorQTL, and performing fine-mapping with SuSiE.

All workflows referenced here are implemented as WDLs in **WARP**.

The *original versions* of these workflows were either created by the GTEx Consortium (see their [GTEx GitHub repository](https://github.com/broadinstitute/gtex-pipeline/tree/master?tab=readme-ov-file)) or the lab for  **Dr. Stephen Montgomery Lab** at Stanford University, with major contributions from **Evin Padhi** and **Jon Nguyen**. Their work formed the foundation for the integrated analysis pipeline described here. Portions of the logic originated from the publicly available repository:

* **AoU-Multiomics-Analysis**
  [https://github.com/AoU-Multiomics-Analysis](https://github.com/AoU-Multiomics-Analysis)

This README explains how the pieces fit together and what each WDL component produces.

---

# Table of Contents

1. [Overview](#overview)
2. [Input Requirements](#input-requirements)
3. [Analysis Flow](#analysis-flow)

   * [1. Ancestry Grouping & Sample Lists](#1-ancestry-grouping--sample-lists)
   * [2. Genotype Preparation (`Prepare_VCF`)](#2-genotype-preparation-prepare_vcf)
   * [3. Genotype Dosage Calculation](#3-genotype-dosage-calculation)
   * [4. RNA Alignment, Counts and Splicing BED](#4-rna-alignment-counts-and-splicing-bed)
   * [5. RNA Phenotype Preparation (`Prepare_eQTL`)](#5-rna-phenotype-preparation-prepare_eqtl)
   * [6. Covariate Creation (`MergeCovariates`)](#6-covariate-creation-mergecovariates)
   * [7. cis-eQTL Mapping (TensorQTL)](#7-cis-eqtl-mapping-tensorqtl)
   * [8. FDR Recalculation & Fine-Mapping Prep](#8-fdr-recalculation--fine-mapping-prep)
   * [9. SuSiE Fine-Mapping (`SusieR`)](#9-susie-fine-mapping-susier)
   * [10. Allele Frequency Calculation](#10-allele-frequency-calculation)
   * [11. SuSiE Aggregation](#11-susie-aggregation)
4. [sQTL Workflow](#sqtl-workflow)
5. [Acknowledgements](#acknowledgements)

---

# Overview

This pipeline generates all inputs, outputs, and intermediate metadata required for a full cis-expression QTL(eQTL) and cis-splicing QTL (sQTL) analysis:

* Genotype preprocessing and pruning
* PLINK files and genotype principal components
* Dosage matrices per ancestry
* Expression and splicing phenotype matrices
* Phenotype PCs and additional grouping metadata
* Covariate tables
* TensorQTL cis-QTL results
* SuSiE fine-mapping outputs
* Aggregated credible sets

---

# **Input Requirements**

To run the workflows described here, you need:

* A joint-called **VCF** containing the relevant samples
* **Research IDs** partitioned by ancestry or subpopulation
* RNA expression quantifications (for eQTLs)
* BAM/CRAM files for splice junction extraction (for sQTLs)
* Associated metadata (sample-level phenotype table)

---

# **Analysis Flow**

The sections below describe each workflow, its purpose, and expected outputs.

---

## 1. Ancestry Grouping & Sample Lists

Prepare a table listing sample IDs for each ancestry/subpopulation.

Outputs:

* Sample lists per group
* Updated tables of sample metadata
* Input tables needed for downstream WDLs

This step is required before running genotype or phenotype workflows per ancestry.

---

## 2. Genotype Preparation (`Prepare_VCF`)

The [Prepare_VCF](https://dockstore.org/workflows/github.com/AoU-Multiomics-Analysis/prepare_QTL/prepare_VCF:develop?tab=info) WDL performs:

* Variant pruning
* Conversion of the VCF to PLINK (`pgen`, `psam`, `pvar`)
* Computation of **genotype PCs**

Outputs:

* Pruned VCF
* PLINK genotype files
* Genotype principal component matrix

These outputs are used for both eQTL and sQTL pipelines.

---

## 3. Genotype Dosage Calculation

The [CalculateGenotypeDosage](https://github.com/AoU-Multiomics-Analysis/prepare_QTL/blob/main/workflows/calculateGenotypeDosage.wdl) WDL generates genotype dosages per ancestry group.

Outputs:

* Two dosage files per ancestry (variant-by-sample dosage matrices)

Because this step uses only the VCF, it may be integrated with `Prepare_VCF` in future versions.

---
## 4. RNA Alignment, Counts and Splicing BED
The [rnaseq_aou.wdl](./rnaseq_aou.wdl) was modified from the original GTEx pipeline and run with the GENCODE v48 GTF. The resulting counts were used as input for downstream expression QTL analysis.

For splicing QTL (sQTL) analysis, the resulting duplicated-marked aligned BAMs were used as a input to the [leafcutter_bam_to_juc wdl](./leafcutter_bam_to_junc.wdl).

This created a junction file that was used as input for the [leafcutter_cluster.wdl](./leafcutter_cluster.wdl), which produces a BED file for downstream  sQTL.

## 5. RNA Phenotype Preparation for eQTL (`Prepare_eQTL`)

The [Prepare_eQTL WDL](https://github.com/AoU-Multiomics-Analysis/prepare_QTL/blob/main/workflows/prepare_eQTL.wdl) processes RNA expression data to generate:

* A **BED-format phenotype matrix**
* **Phenotype PCs** for downstream covariate construction

Inputs typically include:

* Expression quantifications (TPM/CPM or counts)
* Sample metadata
* Gene annotations

Outputs are formatted to match TensorQTL requirements.

---

## 6. Covariate Creation (`MergeCovariates`)

The [MergeCovariates WDL](https://github.com/AoU-Multiomics-Analysis/prepare_QTL/blob/main/workflows/MergeCovariates.wdl) merges:

* Genotype PCs
* Phenotype PCs (expression or splicing)
* Optional grouping variables (for sQTLs)

Outputs:

* A single covariates file for use in TensorQTL

This step ensures consistent ordering and formatting across all inputs.

---

## 7. cis-eQTL Mapping (TensorQTL)

The [TensorQTL cis permutations WDL](https://dockstore.org/workflows/github.com/AoU-Multiomics-Analysis/tensorQTL_cis_permutations:main?tab=info) runs **TensorQTL cis-permutation** mode to compute:

* Nominal associations
* Permutation-based cis-eQTL statistics
* Beta values, empirical p-values, and effect directions

Notes:

* The optional phenotype groups file is **not** required for eQTL analysis
* Results can be written directly into a structured output directory or table
* These outputs form the basis for fine-mapping

---

## 8. FDR Recalculation & Fine-Mapping Prep

After TensorQTL completes:

* Recalculate FDR
* Filter to **FDR ≤ 0.05**
* Format results into a SuSiE-ready table

This step typically involves:

* Computing q-values
* Generating a list of significant gene–variant pairs
* Preparing SuSiE input metadata, including:

  * Expression ID
  * Genomic window coordinates
  * Output prefix names (must match SuSiE input requirements)

---

## 9. SuSiE Fine-Mapping (`SusieR`)

This [SusieR WDL](./susieR_workflow.wdl) performs SuSiE fine-mapping for each cis-window.

Inputs:

* Dosage matrices (from Step 3)
* Significant TensorQTL hits (from Step 7)
* Expression or splicing phenotype metadata
* A consistent `OutputPrefix` per phenotype

Outputs include:

* SuSiE credible sets
* Variant posterior inclusion probabilities
* Fine-mapped credible intervals

Tips:

* Preemptible VMs can be used to reduce cost
* For reproducibility, a pinned Docker SHA is recommended

---

## 10. Allele Frequency Calculation

The [CalculateAF](https://dockstore.org/workflows/github.com/AoU-Multiomics-Analysis/prepare_QTL/calculateAF:main?tab=info) WDL calculates allele frequencies using PLINK.

Outputs:

* Per-variant allele frequencies
* Additional variant summary metrics

This step is optional but useful for interpretation and downstream reporting.

---

## 11. SuSiE Aggregation

The [AggregateSusie WDL](./AggregateSusieWorkflow.wdl) aggregates fine-mapping results across all phenotypes.

Inputs:

* Paths to all SuSiE parquet outputs

  * Use the **“SusieParquet”** (fine-mapped) files
  * Do **not** use the "Full" parquets (contain all tested variants)

Outputs:

* Combined table of all credible sets
* Aggregated fine-mapping metadata
* Summary tables for downstream QTL interpretation

---

# **sQTL Workflow**

The sQTL pipeline shares genotype components with the eQTL workflow but differs in phenotype preparation and covariate structure.

---

## **1. Leafcutter Junc and Cluster Generation**

Run:

* **Bam2Junc** to extract junctions
* **Cluster** to identify splice clusters

Outputs:

* Junc files
* Cluster definitions
* Leafcutter BED files (cluster-level)

---

## **2. Prepare sQTL Phenotypes (`prepare_sQTL`)**

This WDL:

* Consumes Leafcutter BED files
* Generates a **splicing phenotype BED**
* Computes **phenotype PCs**

Older versions of the preprocessing script also produced:

* **PhenotypeGroups** (required for TensorQTL)

This is included as a separate workflow below.

---

## **3. Calculate Phenotype Groups**

If phenotype groups are not emitted by the updated splicing phenotype WDL, a supplementary WDL can generate them.

Outputs:

* PhenotypeGroups file

---

## **4. Merge Covariates (sQTL)**

Identical to the eQTL covariate merging step, but includes:

* Genotype PCs
* Splicing phenotype PCs
* PhenotypeGroups file

Outputs:

* Covariate file for TensorQTL sQTL analysis

---

## **5. TensorQTL cis-sQTL**

This step uses:

* PLINK genotype files
* Splicing BED phenotype matrix
* Covariates
* PhenotypeGroups

Outputs:

* cis-sQTL nominal and permutation results
* Per-cluster association statistics

Downstream fine-mapping can be performed using the same SuSiE workflow if desired.

---

# **Acknowledgements**

This pipeline builds upon extensive work by the **Stephen Montgomery Lab** at Stanford University.
Special thanks to:

* **Evin Padhi**
* **Jon Nguyen**

for developing foundational versions of many scripts and workflows used in this analysis.

Additional integration, optimization, and workflow migration were performed by the All of Us DRC Multiomics and Pipeline Development teams as part of the WARP workflow suite.

