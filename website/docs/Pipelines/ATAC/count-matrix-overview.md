---
sidebar_position: 2
---

# ATAC Count Matrix Overview

The ATAC pipeline's default count matrix output is a h5ad file generated using [SnapATAC2](https://github.com/kaizhang/SnapATAC2) and [AnnData](https://anndata.readthedocs.io/en/latest/index.html). 

<!--- paragraph about contents of matrix --->

<!--- paragraph about anndata.uns, anndata.obs, and .var --->


## Table 1. Global attributes

The global attributes (unstuctured metadata) in the h5ad apply to the whole file, not any specific part. 

| Attribute | Details |
| --- | --- |
| `reference_sequences` | Data frame containing the chromosome sizes for the genome build (i.e., hg38); created using the `chrom_sizes` pipeline input. |


## Table 2. Cell metrics

| Cell Metrics | Program | Details |
|---|---|--------------------|
| `tsse` | [SnapATAC2](https://github.com/kaizhang/SnapATAC2) | Transcription start site enrichment (TSSe) score; lower scores suggest poor data quality. Read more about TSSe scores in the [Definitions section](#definitions) below. |
| `n_fragment` | [SnapATAC2](https://github.com/kaizhang/SnapATAC2) |
| `frac_dup` | [SnapATAC2](https://github.com/kaizhang/SnapATAC2) |
| `frac_mito` | [SnapATAC2](https://github.com/kaizhang/SnapATAC2) |
| `gex_barcodes` | [AnnData](https://anndata.readthedocs.io/en/latest/index.html) |


## Definitions
