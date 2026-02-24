# Admixture & Ancestry Estimation Workflows (WDL)

This directory contains three WDL workflows used for preparing genotype data and estimating genetic ancestry proportions. These workflows are designed to operate on genotype data derived from a **combined reference panel of 1000 Genomes Project (1KG) and Human Genome Diversity Project (HGDP)** samples.

## Background & Interpretation Notes

The 1KG + HGDP reference panel provides broad global coverage but has **known limitations**, including uneven population representation and limited resolution for certain ancestries. As a result:

* Some ancestry components may be **overestimated or underestimated**, depending on the method and reference composition.
* Outputs from these workflows should be **interpreted cautiously** and in the context of these reference limitations.
* Results are best suited for **population-level summaries** rather than precise individual-level ancestry inference.

To address specific biases observed with Rye-based admixture estimation, an alternative ADMIXTURE-based workflow is included and was used to generate a **pruned reference for All of Us (AoU) admixture analyses**.

---

## Workflow 1: `convert_vcf_to_plink_bed`

### Description

Converts a merged VCF file into PLINK binary format (`.bed/.bim/.fam`) for downstream ancestry and admixture analyses.

### Required Inputs

* `prefix` (String): Prefix for output files
* `merged_vcf_shards` (File): Merged VCF file
* `merged_vcf_shards_idx` (File): Index file for the merged VCF

### Tasks & Software

* **Task:** `convert_vcf_to_plink_bed`
* **Software:**

  * PLINK (`/app/bin/plink`)
  * Docker image: `mussmann/admixpipe:3.0`

PLINK is run with:

* `--vcf` input
* `--make-bed` to generate binary PLINK files
* `--double-id` and `--allow-extra-chr` for compatibility with reference data

### Outputs

* `*.bed`: PLINK binary genotype file
* `*.bim`: Variant information file
* `*.fam`: Sample metadata file

---

## Workflow 2: `run_admixture_est_rye`

### Description

Generates ancestry (admixture) estimates using the **Rye** tool, based on PCA-derived eigenvalues and eigenvectors. Outputs are analogous to ADMIXTURE `Q` and `fam` files.

This workflow is useful for fast, PCA-informed ancestry estimation but was observed to **overestimate certain populations** when using the AoU reference.

### Required Inputs

* `eigenvalues_file` (File): PCA eigenvalues
* `eigenvec_file` (File): PCA eigenvectors
* `pop2group_file` (File): Population-to-group mapping
* `prefix` (String): Output file prefix
* `pcs` (Int, default = 20): Number of principal components
* `rounds` (Int, default = 200): Optimization rounds
* `iter` (Int, default = 100): Optimization iterations
* `cpus` (Int, default = 16)
* `docker_image` (String): Rye Docker image

### Tasks & Software

* **Task:** `run_rye`
* **Software:**

  * Rye (`rye.R`)
  * Docker image:
    `us-central1-docker.pkg.dev/broad-dsde-methods/aou-auxiliary/rye-admixture-estimation-tool:v1.0`

### Outputs

* `*.Q`: Admixture proportion file
* `*.fam`: Sample metadata file

Output files are renamed for consistency:

```
<prefix>-<pcs>.Q
<prefix>-<pcs>.fam
```

---

## Workflow 3: `run_admixture`

### Description

Runs the **ADMIXTURE** software directly on PLINK binary files to estimate ancestry proportions.

This workflow was used to generate a **pruned AoU admixture reference**, specifically because **Rye-based estimates were found to overestimate certain populations**, while ADMIXTURE tends to **underestimate those same groups**, providing a complementary and more conservative estimate.

### Required Inputs

* `bed` (File): PLINK `.bed` file
* `bim` (File): PLINK `.bim` file
* `fam` (File): PLINK `.fam` file
* `K_in` (Int, optional): Number of ancestry components (default = 6)
* `num_cpus_in` (Int, optional): Number of CPU threads (default = 4)
* `mem_gb` (Int, default = 120): Memory allocation

### Tasks & Software

* **Task:** `run_admixture`
* **Software:**

  * ADMIXTURE (`/app/bin/admixture`)
  * Docker image: `mussmann/admixpipe:3.0`

ADMIXTURE is executed with multithreading (`-j`) for performance.

### Outputs

* `*.Q`: Ancestry proportion matrix
* `*.P`: Population allele frequency matrix

Output naming follows ADMIXTURE conventions:

```
<basename>.<K>.Q
<basename>.<K>.P
```

---

## Summary

Together, these workflows support:

1. Conversion of VCF data to PLINK format
2. PCA-based admixture estimation using Rye
3. Direct admixture estimation using ADMIXTURE for reference refinement

When interpreting results, users should consider:

* Reference panel composition (1KG + HGDP)
* Method-specific biases (Rye vs. ADMIXTURE)
* The intended use case (population-level inference vs. individual ancestry)

