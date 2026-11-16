# All of Us Phased Small Callsets Pipeline

Builds the three small callsets (Exome, ClinVar, ACAF threshold) from **phased, per-chromosome VCFs**
instead of a VDS, and exports **per-chromosome VCFs** as the priority-one deliverable.

This is a sibling of [`../small_callsets`](../small_callsets), which builds the same callsets from a VDS.
For the callset definitions, BED-construction notes and the high-level analysis overview, see the
WARP docs site:
[All of Us / Small Callsets overview](https://broadinstitute.github.io/warp/docs/All_of_Us/Small_Callsets/overview).

## Run order

| Step | Workflow | Input | Output |
| --- | --- | --- | --- |
| 1 | `sc1_create_basis_mt_from_vcf` | Phased per-chromosome VCF + 3 UCSC BED files | Basis MatrixTable annotated with BED membership, plus a JSON manifest |
| 2 | `sc2_export_basis_mt_to_dense_bed_mts` | SC1 manifest | One dense MatrixTable per BED |
| 3 | `sc4_export_dense_bed_MTs_to_split_dense_MTs` | SC2 manifest | One split (biallelic) MatrixTable per BED |
| 4 | `sc7_export_split_dense_bed_MTs_to_vcf` | SC4 manifest | One VCF per BED, with a phased-VCF header |
| 5 | `CreateVcfIndex` ([`../phasing`](../phasing)) | SC7 VCFs | Tabix indexes |

There is no SC3, SC5 or SC6 here. SC5 (BGEN) and SC6 (PLINK) are deliberately excluded: neither format
can represent phase, so both would silently discard it.

## Relationship to `../small_callsets`

SC1-from-VCF and SC7 are new and exist only here. **SC2 and SC4 are forked copies** of their VDS-pipeline
siblings — byte-identical as of 2026-11-16 — so the phased path can change without touching VDS
production. Each forked file carries a `PHASED FORK` header comment. **Fix bugs in both copies** unless
the fix is genuinely phasing-specific.

Because Dockstore workflow names must be unique across the repo, the forked workflows are registered as
`sc2_export_basis_mt_to_dense_bed_mts_phased` and `sc4_export_dense_bed_MTs_to_split_dense_MTs_phased`,
even though the file names match their siblings.

## Changelogs

- [`sc1_create_basis_mt_from_vcf.changelog.md`](sc1_create_basis_mt_from_vcf.changelog.md)
- [`sc2_export_basis_mt_to_dense_bed_mts.changelog.md`](sc2_export_basis_mt_to_dense_bed_mts.changelog.md)
- [`sc4_export_dense_bed_MTs_to_split_dense_MTs.changelog.md`](sc4_export_dense_bed_MTs_to_split_dense_MTs.changelog.md)
