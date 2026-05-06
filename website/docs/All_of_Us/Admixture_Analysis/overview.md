---
title: Admixture & Ancestry Estimation Workflows (WDL)
className: aou-doc-page
---

<div className="aou-folder-text">

This section describes the All of Us admixture workflows used to estimate global ancestry proportions from the ancestry pipeline outputs.

## Quick Summary

* **Purpose:** Support two analysis paths (Admixture Rye and Admixture Unsupervised) for admixture estimation in All of Us.
* **Admixture Rye:** `run_preprocess_admixture_est_rye` → `run_admixture_est_rye`.
* **Admixture unsupervised:** `convert_vcf_to_plink_bed` → `run_admixture` (unsupervised clustering).

## Background and Interpretation Notes

This pipeline uses the Gnomad 3.1.2 reference panel (1KG + HGDP). This panel provides broad global coverage but has known limitations, including uneven population representation and limited resolution for some ancestries. As a result:

* Some ancestry components may be overestimated or underestimated, depending on method and reference composition.
* Outputs should be interpreted cautiously and in context of reference limitations.
* Results are best suited for population-level summaries rather than precise individual-level ancestry inference.

---

## Admixture Rye Analysis

For the All of Us Admixture Rye analysis, run these workflows in order:

1. `run_preprocess_admixture_est_rye`
2. `run_admixture_est_rye`

### Workflow 1: `run_preprocess_admixture_est_rye`

#### Description

Generates Rye-ready input files from ancestry pipeline outputs. This workflow prepares:

* Rye eigenvalues file
* Rye eigenvector file
* Rye population-to-group mapping file

It is designed to consume PCA/ancestry outputs generated upstream in the ancestry workflow.

#### Required Inputs

* `eigenvalues_url` (String): Path to PCA eigenvalues
* `training_pca_url` (String): Path to training PCA table
* `ancestry_data_url` (String): Path to ancestry prediction output containing projected PCA features
* `prefix` (String): Output prefix
* `cpus` (Int, default = 16)
* `docker_image` (String, default = `hailgenetics/hail:0.2.67`)

#### Outputs

* `rye_eigenval` (`<prefix>_rye.eigenvalues`)
* `rye_eigenvec` (`<prefix>_rye.eigenvec`)
* `rye_pop2group` (`<prefix>_rye.pop2group`)

### Workflow 2: `run_admixture_est_rye`

#### Description

Runs the Rye tool to estimate ancestry proportions using PCA-derived eigenvalues/eigenvectors and population grouping metadata from preprocessing.

#### Required Inputs

* `eigenvalues_file` (File): Rye-compatible eigenvalues file
* `eigenvec_file` (File): Rye-compatible eigenvectors file
* `pop2group_file` (File): Population-to-group mapping file
* `prefix` (String): Output file prefix
* `pcs` (Int, default = 20)
* `rounds` (Int, default = 200)
* `iter` (Int, default = 100)
* `cpus` (Int, default = 16)
* `docker_image` (String): Rye Docker image

#### Outputs

* `*.Q`: Admixture proportion file
* `*.fam`: Sample metadata file

Output files are renamed for consistency:

```
<prefix>-<pcs>.Q
<prefix>-<pcs>.fam
```

---

## Admixture Unsupervised (for underrepresented populations)

Running Admixture unsupervised allow for further customization of the reference panel that may better account for underrepresented populations. Run these workflows in order:

1. `convert_vcf_to_plink_bed`
2. `run_admixture` (unsupervised clustering)

### Workflow 3: `convert_vcf_to_plink_bed`

#### Description

Converts merged VCF training inputs from the ancestry pipeline into PLINK binary format (`.bed/.bim/.fam`) for downstream ADMIXTURE clustering.

#### Required Inputs

* `prefix` (String): Prefix for output files
* `merged_vcf_shards` (File): Merged VCF file
* `merged_vcf_shards_idx` (File): Index file for the merged VCF

#### Tasks & Software

* **Task:** `convert_vcf_to_plink_bed`
* **Software:**
  * PLINK (`/app/bin/plink`)
  * Docker image: `mussmann/admixpipe:3.0`

PLINK is run with:

* `--vcf` input
* `--make-bed`
* `--double-id` and `--allow-extra-chr`

#### Outputs

* `*.bed`: PLINK binary genotype file
* `*.bim`: Variant information file
* `*.fam`: Sample metadata file

### Workflow 4: `run_admixture`

#### Description

Runs ADMIXTURE in unsupervised mode on PLINK binary files to estimate ancestry components and support downstream pruned-reference construction.

#### Required Inputs

* `bed` (File): PLINK `.bed` file
* `bim` (File): PLINK `.bim` file
* `fam` (File): PLINK `.fam` file
* `K_in` (Int, optional, default = 6): Number of ancestry components
* `num_cpus_in` (Int, optional, default = 4): Number of CPU threads
* `mem_gb` (Int, default = 120): Memory allocation

#### Tasks & Software

* **Task:** `run_admixture`
* **Software:**
  * ADMIXTURE (`/app/bin/admixture`)
  * Docker image: `mussmann/admixpipe:3.0`

#### Outputs

* `*.Q`: Ancestry proportion matrix
* `*.P`: Population allele frequency matrix

Output naming follows ADMIXTURE conventions:

```
<basename>.<K>.Q
<basename>.<K>.P
```

---

## Notes for Interpretation

When interpreting outputs, consider:

* Reference panel composition (1KG + HGDP)
* Method-specific behavior (Rye vs. unsupervised ADMIXTURE)
* Intended use case (population-level inference vs. individual-level calls)

</div>

