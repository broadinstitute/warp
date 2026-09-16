---
sidebar_position: 6
slug: /Pipelines/Optimus_Pipeline/v4_chemistry
---

# 10x 3' v4 (GEM-X) Chemistry Support

This page explains how to configure and run Optimus with 10x Genomics 3' v4 (GEM-X) chemistry, covering FASTQ inputs, reference files, barcode whitelist selection, and counting mode.

## Overview

Optimus natively supports 10x 3' v4 (GEM-X) chemistry by setting `tenx_chemistry_version = 4`. The v4 chemistry uses the same cell barcode and UMI lengths as v3:

| Feature | v2 | v3 | v4 (GEM-X) |
| --- | --- | --- | --- |
| Cell barcode length | 16 bp | 16 bp | 16 bp |
| UMI length | 10 bp | 12 bp | 12 bp |
| R1 total length | 26 bp | 28 bp | 28 bp |
| Whitelist (GCP) | `737K-august-2016.txt` | `3M-february-2018.txt` | `3M-3pgex-may-2023.txt` or `3M-3pgex-may-2023_TRU.txt` |

## FASTQ Inputs

Three FASTQ inputs are required for v4 chemistry. All lanes must be provided as arrays:

| Input | Contents | Required for v4? |
| --- | --- | --- |
| `Optimus.r1_fastq` | R1 reads: cell barcode + UMI | Yes |
| `Optimus.r2_fastq` | R2 reads: cDNA sequence | Yes |
| `Optimus.i1_fastq` | I1 index reads: sample barcodes | **Yes — pipeline will fail if omitted** |

## Reference Files

The STAR reference tarball and GTF must be from matched builds. For a full list of available GCP-hosted references (human GRCh38, mouse GRCm39, and others), see the [BuildIndices example references](https://broadinstitute.github.io/warp/docs/Pipelines/BuildIndices_Pipeline/README#example-references).

:::warning Match your reference and GTF
Always use a STAR index and GTF from the same genome build and annotation release. Mismatched files will produce incorrect gene assignments.
:::

:::note Pseudogene handling
WARP GENCODE-based references are not filtered to remove pseudogenes, unlike 10x Genomics Cell Ranger references. This can result in small count differences for multi-mapped pseudogene reads. See the [pseudogene handling note](./README.md#pseudogene-handling) for details.
:::

## Chemistry and Counting Mode

Set `tenx_chemistry_version = 4` for all 10x v4 (GEM-X) data.

For single-nucleus RNA-seq (the most common use case with v4), also set:

```json
"Optimus.counting_mode": "sn_rna"
```

`sn_rna` mode uses the `GeneFull_Ex50pAS` STARsolo feature, which counts reads aligned to whole transcripts (introns + exons) to account for nuclear pre-mRNA. For single-cell (not nucleus) v4 data, use `"sc_rna"` instead.

## Barcode Whitelist Selection

:::caution This choice affects cell calling
The whitelist you select directly determines which barcodes are recognized as valid cells and affects comparability with Cell Ranger outputs. Use the correct version for your kit.
:::

Because 10x Genomics updated the v4 barcode list between Cell Ranger releases, Optimus exposes a `tenx_chemistry_subversion` input:

| `tenx_chemistry_subversion` | Whitelist GCS path | Cell Ranger version | When to use |
| --- | --- | --- | --- |
| `"v4_TRU"` *(default when unspecified)* | `gs://gcp-public-data--broad-references/optimus_whitelists/3M-3pgex-may-2023_TRU.txt` | v9.0 and later | **Recommended for all new analyses** |
| `"v4"` | `gs://gcp-public-data--broad-references/optimus_whitelists/3M-3pgex-may-2023.txt` | v8.0 and v8.0.1 | Only when backward compatibility with CR v8 outputs is required |

For a complete list of all available whitelists, see the [Barcode whitelist options](./README.md#barcode-whitelist-options) table in the main Optimus documentation.

## Example Inputs

Downsampled plumbing test examples (human GRCh38, two lanes) are available in the repository:
- [human_v4_TRU_example.json](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/optimus/example_inputs/human_v4_TRU_example.json) — Cell Ranger v9.0+ (recommended)
- [human_v4_example.json](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/optimus/example_inputs/human_v4_example.json) — Cell Ranger v8.0/v8.0.1

For single-nucleus data, add `"Optimus.counting_mode": "sn_rna"` to your input JSON. For Cell Ranger v8.0/v8.0.1 backward compatibility, change `tenx_chemistry_subversion` to `"v4"`.

## Read Structure

STARsolo is run with:

```
--soloCBstart 1 --soloCBlen 16
--soloUMIstart 17 --soloUMIlen 12
```

R1 is expected to be 28 bp (16 bp CB + 12 bp UMI). The pipeline validates this length by default; set `ignore_r1_read_length = true` to skip validation if your R1 length differs.

## Validation

After running, verify that:
- The `whitelist_input_used` output reflects the expected `3M-3pgex-may-2023*.txt` path.
- STARsolo Summary.csv shows CB length 16 and UMI length 12.
- Per-gene and per-cell Spearman correlations with Cell Ranger outputs are ≥ 0.97 (expected with matched references).
