### Pipelines for Generating Smaller Callsets from a Joint WGS VDS

This directory contains a staged small-callset pipeline that starts from a joint WGS VDS and produces bed-restricted deliverables (exome, ClinVar, ACAF threshold), including Hail MT outputs and downstream BGEN / PLINK exports.

Most stages are implemented as a **WDL wrapper** that launches a Dataproc/Hail job and executes a paired **Python script** for data transformation logic.

## Workflow map: WDL + Python combinations

| Stage | WDL | Python Script | High-level purpose | Main output / handoff |
|---|---|---|---|---|
| SC1 | `sc1_create_basis_mt_from_vds.wdl` | `sc1_create_basis_mt_from_vds.py` | Read joint VDS, filter to union of target UCSC BED intervals, create a basis dense MT with bed membership flags. | Basis MT JSON manifest (`basis_mt_json_file`) pointing to basis MT path. |
| SC2 | `sc2_export_basis_mt_to_dense_bed_mts.wdl` | `sc2_export_basis_mt_to_dense_bed_mts.py` | Split basis MT into one dense MT per bed category (`exome`, `clinvar`, `acaf_threshold`). | Dense-MT JSON manifest (`dense_mts_json_file`). |
| SC4 | `sc4_export_dense_bed_MTs_to_split_dense_MTs.wdl` | `sc4_export_dense_bed_MTs_to_split_dense_MTs.py` | Convert dense multi-allelic MTs into split (bi-allelic) dense MTs and recompute variant summary fields. | Split-dense-MT JSON manifest (`split_dense_mts_json_file`). |
| SC5 | `sc5_export_split_dense_bed_MTs_to_bgen.wdl` | `sc5_export_split_dense_bed_MTs_to_bgen.py` | Export split dense MTs to chromosome-level BGEN outputs per bed category. | BGEN/sample manifests (e.g., `exome_manifest`) and BGEN files in GCS. |
| SC5.2 | `sc5_2_bgenix_index.wdl` | _(none)_ | Index exported BGEN files using `bgenix` (`.bgi`). | `.bgi` files and index FOFN. |
| SC6 | `sc6_export_split_dense_bed_MTs_to_plink_bed.wdl` | `sc6_export_split_dense_bed_MTs_to_plink_bed.py` | Export split dense MTs to chromosome-level PLINK BED deliverables per bed category. | PLINK manifest(s) (e.g., `.bim` manifest) and PLINK files in GCS. |

## Recommended run order and dependencies

1. Run **SC1** to generate the basis MT manifest JSON.
2. Run **SC2** using SC1’s JSON output.
3. Run **SC4** using SC2’s JSON output.
4. Run **SC5** and/or **SC6** using SC4’s JSON output, depending on required deliverable format.
5. Run **SC5.2** after SC5 to build BGEN indexes.

In short:  
`SC1 -> SC2 -> SC4 -> (SC5 -> SC5.2) and/or SC6`

## How BED files were created

### ACAF BED (All Common Alleles/Frequencies)

This BED is used to subset the full callset (500K samples) to variants that are common in at least one ancestry group or common overall for ACAF-level analyses.

**Inclusion criteria (from VAT in BigQuery):**
- A variant is included if **either** of these is true:
  - Per-ancestry allele count threshold: `gvs_{ancestry}_ac >= 100` for **any** ancestry (OR across ancestries)
  - Overall max AF threshold: `gvs_max_af > 0.01`

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

- BEDs are position-based (`contig`, `position`) and do not preserve allele-level specificity.
- At multi-allelic sites, if one allele qualifies, interval-based downstream filtering keeps the full site position.
- The `AC >= 100` threshold is an absolute allele count (not frequency), so effective frequency threshold varies by ancestry sample size.
- `bedtools merge` here uses exact overlap/adjacency behavior (no extra gap bridging without `-d`).

## Notes

- “BED” in this pipeline primarily refers to **UCSC BED interval files** used for region filtering.
- PLINK BED output in SC6 is a **genotype file format** and is different from UCSC BED intervals.
- Several WDLs are designed for Terra/Dataproc execution and pass Spark/Hail settings through to the paired Python scripts.