# sc1_create_basis_mt_from_vcf.py -- validation notes

Interactive validation performed on the hail fat box against real phased AoU data, before wiring this
script into a WDL. Recorded here so the testing doesn't only live in a Jupyter session.

## Test region

Liz identified `chr21:10494999-10547000` as a region common to all 3 BED files (exome, clinvar,
acaf_threshold), which makes it a good single window for exercising the whole BED-annotate/filter step
at once instead of guessing at coordinates.

Input: `gs://prod-drc-broad/aou_phasing/folf/chr21.aou.v9.phased.filtered.vcf.bgz`

BED files used (note: `exome_bed_path` and `acaf_threshold_bed_path` were initially swapped by mistake --
confirmed by inspecting each file's interval sizes: the `ac_af_cutoff_*` file has single-variant-sized
1-3bp intervals, consistent with an AC/AF-cutoff-selected site list, while the `gencode_*` file has
exon-sized intervals of tens to thousands of bp, consistent with an exome capture/gene-model BED):
* `exome_bed_path`: `gencode.v42.primary_assembly.basic.annotation.gtf.bedops.cut.slop15.sorted.merged.pgx.chrM_last.bed`
* `clinvar_bed_path`: `clinvar_foxtrot_v4.sorted.ucsc.bed`
* `acaf_threshold_bed_path`: `ac_af_cutoff_foxtrot_v4.sorted.ucsc.bed`

## 1. BED annotate/filter runs against all 3 BEDs and finds real overlap

```python
mt.aggregate_rows(hl.struct(
    exome=hl.agg.count_where(mt.exome),
    clinvar=hl.agg.count_where(mt.clinvar),
    acaf_threshold=hl.agg.count_where(mt.acaf_threshold),
))
```
Confirmed non-zero hits for all 3 BED types in this region after correcting the path swap above (an
earlier run, before the swap was caught, returned `Struct(exome=1134, clinvar=2, acaf_threshold=379)` --
those specific counts were against the swapped/mislabeled paths and should not be trusted as final
numbers; rerun this check and record the corrected counts here once available).

## 2. Phase survives import + variant_qc/info recompute (make_basis_mt_from_vcf)

```python
entries = mt.entries()
het_entries = entries.filter(entries.GT.is_het())
het_entries.select(gt=het_entries.GT, is_phased=het_entries.GT.phased).show(5)
```
Result: real heterozygous calls came back phased (`0|1`/`1|0`, `is_phased=True`) after passing through
`hl.import_vcf` -> BED annotate/filter -> `hl.variant_qc` -> the `info` struct merge in
`make_basis_mt_from_vcf`. This is the core guarantee this script exists to provide.

## 3. split_multi_hts phase preservation (downstream SC4 step, tested here since this is where the
multiallelic site was found)

Test site: `chr21:10495314`, `alleles=["G","A","C","T"]` (info.AC = [214, 1, 1] -- A is common, C and T
are each carried by exactly one person in this cohort).

**3a. Two different single-alt-vs-ref carriers, opposite haplotype sides:**
* `1686377`, original GT `0|3` (T on the right haplotype)
* `2647058`, original GT `2|0` (C on the left haplotype)

```python
mt_locus = mt.filter_rows(mt.locus == hl.parse_locus("chr21:10495314", reference_genome="GRCh38"))
mt_split = hl.split_multi_hts(mt_locus)
split_entries = mt_split.entries()
check_samples = ["1686377", "2647058", "1009368"]  # + a common A-carrier for contrast
filtered = split_entries.filter(hl.literal(set(check_samples)).contains(split_entries.s))
filtered.select(a_index=filtered.a_index, gt=filtered.GT, is_phased=filtered.GT.phased).show(9)
```

| sample | orig GT | a_index=1 (A) | a_index=2 (C) | a_index=3 (T) |
|---|---|---|---|---|
| 1009368 | `1\|0` | `1\|0` | `0\|0` | `0\|0` |
| 1686377 | `0\|3` | `0\|0` | `0\|0` | `0\|1` (T stayed right) |
| 2647058 | `2\|0` | `0\|0` | `1\|0` (C stayed left) | `0\|0` |

All rows `is_phased=True`. Each carrier's real alt landed in the correct split row, on the same
haplotype side it started on; nobody else's allele leaked into a split row they don't belong in.

**3b. a_index-to-allele mapping is fixed and consistent** (not something that varies row to row):
```python
rows = mt_split.rows()
rows.select(a_index=rows.a_index).show()
```
Result: `a_index=1` -> `["G","A"]`, `a_index=2` -> `["G","C"]`, `a_index=3` -> `["G","T"]`, always.

**3c. Synthetic case: one person carrying two different alts at once** (not naturally present at this
site, since C and T are each carried by only one person total -- injected manually to close this gap):
```python
synthetic_sample = "1009368"
mt_synthetic = mt_locus.annotate_entries(
    GT=hl.if_else(mt_locus.s == synthetic_sample, hl.call(1, 2, phased=True), mt_locus.GT)
)
mt_split_synth = hl.split_multi_hts(mt_synthetic)
split_entries = mt_split_synth.entries()
filtered = split_entries.filter(split_entries.s == synthetic_sample)
filtered.select(a_index=filtered.a_index, gt=filtered.GT, is_phased=filtered.GT.phased).show()
```
Synthetic genotype `1|2` (A on the left, C on the right) split correctly: `1|0` / `0|1` / `0|0` for the
A/C/T splits respectively, all `is_phased=True`.

**Conclusion: `split_multi_hts` is safe to use on phased data** for this small callsets pipeline -- all
three at-risk scenarios (simple phased het, two different rare alts on different people, one person with
two different alts) preserve phase and haplotype side correctly.

## Open items not covered by this validation

* `info.HC` was found to duplicate the recomputed `info.homozygote_count` on every row checked --
  original meaning and whether it should be dropped is still an open question for Liz/Jennifer.
* This was tested on one autosomal region only. chrX (haploid non-PAR calls) and higher-ploidy edge
  cases are untested.
* The `AC`/`AF`/`AN`/`homozygote_count` recompute logic itself (not just phase survival) hasn't been
  cross-checked against an independent source for correctness on phased data specifically.

## How to rerun this

See the sc1 module's functions (`import_phased_vcf`, `filter_to_test_intervals_if_needed`,
`annotate_beds_and_filter`, `make_basis_mt_from_vcf`) called directly from a fat-box notebook -- import
the script as a module (`import sc1_create_basis_mt_from_vcf as sc1`) rather than running it as a CLI
script, so intermediate MatrixTables can be inspected between steps.
