---
title: Admixture Estimation Workflows (WDL)
sidebar_position: 1
className: aou-doc-page
---

<div className="aou-folder-text">

This section describes the All of Us admixture workflows used to estimate global ancestry proportions from the ancestry pipeline outputs. There are two supported analysis paths, each consisting of two workflows run in sequence.

## Background

This pipeline uses the GnomAD 3.1.2 reference panel (1KG + HGDP), which provides broad global coverage but has known limitations including uneven population representation and limited resolution for some ancestries. Outputs are best suited for population-level summaries rather than precise individual-level ancestry inference.

---

## Admixture Rye Analysis

The Admixture Rye path uses the [Rye tool](https://github.com/healthdisparities/rye) to estimate ancestry proportions from PCA data. Run these workflows in order:

| Step | Workflow | Description | WDL |
| :--: | --- | --- | :--: |
| 1 | [Admixture Rye Preprocessing](./run_preprocess_admixture_est_rye.md) | Generates Rye-compatible eigenvalues, eigenvectors, and population-to-group mapping files from ancestry pipeline PCA outputs. | [WDL](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/run_preprocess_admixture_est_rye.wdl) |
| 2 | [Admixture Rye](./run_admixture_est_rye.md) | Runs Rye to estimate ancestry proportions; outputs `.Q` and `.fam` files. | [WDL](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/run_admixture_est_rye.wdl) |

---

## Admixture Unsupervised Analysis

The unsupervised path uses [ADMIXTURE](http://dalexander.github.io/admixture/) directly on genotype data, allowing for customization of the reference panel to better account for underrepresented populations. Run these workflows in order:

| Step | Workflow | Description | WDL |
| :--: | --- | --- | :--: |
| 1 | [Admixture Unsupervised Preprocessing](./convert_vcf_to_plink_bed.md) | Converts merged VCF inputs from the ancestry pipeline into PLINK binary format for downstream ADMIXTURE clustering. | [WDL](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/convert_vcf_to_plink_bed.wdl) |
| 2 | [Admixture Unsupervised](./run_admixture.md) | Runs ADMIXTURE in unsupervised mode; outputs ancestry proportion (`.Q`) and allele frequency (`.P`) matrices. | [WDL](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/run_admixture.wdl) |

---

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>

