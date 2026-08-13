---
title: All of Us Small Callsets Pipeline
sidebar_position: 1
className: aou-doc-page
---

<div className="aou-folder-text">

This section documents the All of Us small-callset workflows used to derive smaller, region-restricted callsets from a joint WGS VDS and export downstream analysis-ready formats. In this context, the three primary small callsets are: **Exome** (variants within exome-target BED intervals), **ClinVar** (variants at sites with ClinVar classification annotations), and **ACAF threshold** (common-variant sites selected using AoU VAT-derived thresholds, including per-ancestry `AC >= 100` or overall `max AF > 0.01`, as described in the QC/BED-generation notes).

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

## How BED files were created

### ACAF BED (All Common Alleles/Frequencies)

This BED is used to subset the full callset (500K samples) to variants that are common in at least one ancestry group or common overall for ACAF-level analyses.

**Inclusion criteria (from VAT in BigQuery):**

* A variant is included if **either** of these is true:
  * Per-ancestry allele count threshold: `gvs_{ancestry}_ac >= 100` for **any** ancestry (OR across ancestries)
  * Overall max AF threshold: `gvs_max_af > 0.01`

```python
ac_where_block = "(" + " OR ".join([f'(gvs_{a}_ac >=100)' for a in ancestries]) + ")"

query = f"""
SELECT DISTINCT contig, position
FROM `{bq_table}`
WHERE {ac_where_block}
   OR (gvs_max_af > 0.01)
ORDER BY contig, position
"""
```

### ClinVar BED

ClinVar BED positions are pulled from VAT where `clinvar_classification` is present.

```python
q = f"""SELECT DISTINCT contig, position
FROM `{bq_table}`
WHERE
array_length(clinvar_classification) > 0
ORDER BY contig,position
"""

vat_clinvar_df = (
	bqclient.query(q)
	.result()
	.to_dataframe(
		create_bqstorage_client=True,
		progress_bar_type='tqdm',
	)
)
```

### Shared BED construction steps (ACAF + ClinVar)

Queried positions are converted to 0-based, half-open BED coordinates:

```python
def create_vat_bed_df(vat_df: pd.DataFrame):
	vat_bed_df = vat_df.copy(deep=True)
	vat_bed_df['position_start'] = vat_bed_df['position'] - 1
	vat_bed_df['position_end'] = vat_bed_df['position']
	vat_bed_df = vat_bed_df.drop(columns='position')
	return vat_bed_df
```

Then intervals are sorted and merged:

```python
def sort_and_merge_bed_file(bedfile: str, output_filename: str):
	! echo "Before:" && wc -l {bedfile}
	! ./bedtools sort -i {bedfile} | ./bedtools merge | sort --version-sort > {output_filename}
	! echo "After:" && wc -l {output_filename}
	! cut -f 1 {output_filename} | uniq
```

### Caveats

* BEDs are position-based (`contig`, `position`) and do not preserve allele-level specificity.
* At multi-allelic sites, if one allele qualifies, interval-based downstream filtering keeps the full site position.
* The `AC >= 100` threshold is an absolute allele count (not frequency), so effective frequency threshold varies by ancestry sample size.
* `bedtools merge` here uses exact overlap/adjacency behavior (no extra gap bridging without `-d`).

## Additional Processing Notes

* **BED semantics:** the interval resources used here are position-based and do not preserve allele-level specificity.
* **Multi-allelic behavior:** if one allele at a site qualifies for a BED-derived interval set, downstream interval-based filtering keeps the full site position.
* **ACAF threshold meaning:** `AC >= 100` is an absolute allele-count threshold, not a frequency threshold, so its effective frequency interpretation varies by ancestry sample size.
* **Format naming:** `BED` in the interval-selection context refers to UCSC BED files, while the `PLINK BED` outputs from `SC6` are genotype-format deliverables.

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
