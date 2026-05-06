---
title: Ancestry Analysis Workflows (WDL)
sidebar_position: 1
className: aou-doc-page
---

<div className="aou-folder-text">
This section documents the All of Us ancestry analysis workflows used to convert large callsets, identify shared high-quality sites, infer ancestry labels, and perform downstream cohort QC/relatedness analyses.

## Quick Summary

* **Purpose:** Convert and filter large callsets, build shared high-quality site sets, infer ancestry labels, and perform post-inference sample QC.
* **Primary Outputs:** Per-chromosome VCFs, HQ-intersection callsets, ancestry predictions, relatedness artifacts, and interactive QC plots.

## High-level Analysis Flow

Ancestry analysis is typically run in two phases:

1. **Core ancestry inference preparation and prediction**
2. **Post-inference relatedness and outlier QC**

The first phase produces ancestry labels and PCA features used by the second phase.

## Phase 1: Core Ancestry Inference

Run these workflows in order:

| Step | Workflow | Description | WDL |
| :--: | --- | --- | :--: |
| 1 | [VDS to VCF](./vds_to_vcf.md) | Converts a large VDS into per-contig full and sites-only VCF files for downstream ancestry processing. | [WDL](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/vds_to_vcf.wdl) |
| 2 | [Determine HQ Sites Intersection](./determine_hq_sites_intersection.md) | Intersects training and input-data variants and filters both datasets to the shared high-quality site set. | [WDL](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/determine_hq_sites_intersection.wdl) |
| 3 | [Run Ancestry Inference](./run_ancestry.md) | Trains and applies a PCA-based Random Forest classifier to infer ancestry labels and generate plots. | [WDL](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_ancestry.wdl) |

## Phase 2: Post-inference Cohort QC and Relatedness

Depending on analysis goals, these workflows are used after Phase 1 outputs are available:

| Step | Workflow | Description | WDL |
| :--: | --- | --- | :--: |
| 4 | [Run Relatedness](./run_relatedness) | Computes sample-pair relatedness and a maximal independent set of samples flagged for relatedness deconfounding. | [WDL](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_relatedness.wdl) |
| 5 | [Run Sample Outlier QC](./run_sample_outlier_qc) | Joins ancestry predictions with callset metrics and identifies ancestry-stratified sample QC outliers. | [WDL](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_sample_outlier_qc.wdl) |
| 6 | [Run Sample Outlier QC Plotting](./run_sample_outlier_qc_plotting) | Joins demographics and generates interactive PC and metric visualizations for outlier QC review. | [WDL](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_sample_outlier_qc_plotting.wdl) |

## Notes

* **Data dependency:** `run_sample_outlier_qc` requires ancestry outputs from `run_ancestry`.
* **Plotting dependency:** `run_sample_outlier_qc_plotting` requires outputs from `run_sample_outlier_qc`.
* **Relatedness usage:** `run_relatedness` is commonly used to flag related samples before downstream association analyses.

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>