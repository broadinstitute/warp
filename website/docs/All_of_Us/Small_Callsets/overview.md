---
title: All of Us Small Callsets Pipeline
sidebar_position: 1
className: aou-doc-page
---

<div className="aou-folder-text">

This section documents the All of Us small-callset workflows used to derive smaller, region-restricted callsets from a joint WGS VDS and export downstream analysis-ready formats.

## Quick Summary

* **Purpose:** Generate exome-, ClinVar-, and ACAF-restricted callsets from a joint WGS VDS and export them as Hail MatrixTables, BGEN, and PLINK deliverables.
* **Primary Outputs:** Manifest JSON files for stage-to-stage handoff, per-category dense and split dense MatrixTables, chromosome-level BGEN files with indexes, and chromosome-level PLINK BED outputs.

## Input Requirements

To run the workflows described below, you generally need:

* A joint-called WGS `VDS`
* UCSC `BED` interval resources defining the target regions
* A Dataproc-compatible execution environment for Hail-based processing
* Google Cloud Storage output locations for intermediate and final artifacts
* For BGEN indexing, the list or manifest of exported `.bgen` files from `SC5`

## Ordered Analysis Flow

The table below reflects the typical end-to-end run order for the current Small Callsets analysis path.

| Order | Stage | Workflow / Component | Documentation | WDL | Run next |
| :--: | --- | --- | --- | --- | --- |
| 1 | Basis MT creation | Create basis MatrixTable from VDS | No dedicated page yet | [sc1_create_basis_mt_from_vds.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/small_callsets/sc1_create_basis_mt_from_vds.wdl) | Use the basis MT manifest as input to dense bed MT export. |
| 2 | Dense MT export | Export basis MT to dense bed MatrixTables | No dedicated page yet | [sc2_export_basis_mt_to_dense_bed_mts.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/small_callsets/sc2_export_basis_mt_to_dense_bed_mts.wdl) | Use the dense MT manifest as input to split dense MT export. |
| 3 | Split dense MT export | Export dense bed MatrixTables to split dense MatrixTables | No dedicated page yet | [sc4_export_dense_bed_MTs_to_split_dense_MTs.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/small_callsets/sc4_export_dense_bed_MTs_to_split_dense_MTs.wdl) | Branch to BGEN export and/or PLINK export. |
| 4 | BGEN export | Export split dense MatrixTables to BGEN | No dedicated page yet | [sc5_export_split_dense_bed_MTs_to_bgen.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/small_callsets/sc5_export_split_dense_bed_MTs_to_bgen.wdl) | Run BGEN indexing on exported `.bgen` files. |
| 5 | BGEN indexing | Index BGEN files with `bgenix` | No dedicated page yet | [sc5_2_bgenix_index.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/small_callsets/sc5_2_bgenix_index.wdl) | Consume indexed BGEN outputs downstream. |
| 6 | PLINK export | Export split dense MatrixTables to PLINK BED | No dedicated page yet | [sc6_export_split_dense_bed_MTs_to_plink_bed.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/small_callsets/sc6_export_split_dense_bed_MTs_to_plink_bed.wdl) | Consume PLINK outputs downstream. |

## Practical Run Notes

* **Core path:** `SC1 -> SC2 -> SC4`.
* **BGEN branch:** `SC1 -> SC2 -> SC4 -> SC5 -> SC5.2`.
* **PLINK branch:** `SC1 -> SC2 -> SC4 -> SC6`.
* **Shared dependency:** both export branches depend on the split dense MatrixTable manifests produced by `SC4`.
* **Execution model:** most workflows are WDL orchestration wrappers around Dataproc/Hail execution, with paired Python scripts implementing the transformation logic.

## Additional Processing Notes

* **BED semantics:** the interval resources used here are position-based and do not preserve allele-level specificity.
* **Multi-allelic behavior:** if one allele at a site qualifies for a BED-derived interval set, downstream interval-based filtering keeps the full site position.
* **ACAF threshold meaning:** `AC >= 100` is an absolute allele-count threshold, not a frequency threshold, so its effective frequency interpretation varies by ancestry sample size.
* **Format naming:** `BED` in the interval-selection context refers to UCSC BED files, while the `PLINK BED` outputs from `SC6` are genotype-format deliverables.

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
