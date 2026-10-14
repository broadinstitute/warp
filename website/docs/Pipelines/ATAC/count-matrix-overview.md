---
sidebar_position: 2
---

# ATAC Count Matrix Overview

The [ATAC pipeline's](README.md) default count matrix output is an h5ad file generated using [SnapATAC2](https://github.com/kaizhang/SnapATAC2) and [AnnData](https://anndata.readthedocs.io/en/latest/index.html). 

The h5ad file contains unstructured metadata (`h5ad.uns`; [Table 1](#table-1-global-attributes)) as well as per-barcode quality metrics (`h5ad.obs`; [Table 2](#table-2-cell-metrics)). It also contains an equivalent gene expression barcode for each ATAC barcode. Raw fragments are stored in the `h5ad.obsm['insertion']` property of the h5ad file. For more information, see the [`import_data` function](https://kzhang.org/SnapATAC2/api/_autosummary/snapatac2.pp.import_data.html#snapatac2.pp.import_data) in the SnapATAC2 documentation.

The h5ad file does not contain per-gene metrics, meaning the variables/features data frame (`h5ad.var`) is empty.

The pipeline writes two versions of this file. `snap_metrics` (`<input_id>.metrics.h5ad`) is unfiltered: it has one row for every barcode with at least one fragment. `snap_metrics_filtered` (`<input_id>.filtered.metrics.h5ad`, produced when `filter_cells_to_expected` is true) keeps only the top `atac_expected_cells` barcodes ranked by `n_frag_overlap_peak`. Both carry the `is_cell` and `is_retained` flags described below; in the Multiome pipeline, JoinBarcodes and peak calling use the filtered file.


## Table 1. Global attributes

The global attributes (unstuctured metadata) in the h5ad apply to the whole file, not any specific part. 

| Attribute | Program | Details |
| --- | --- | --- |
| `reference_sequences` | [SnapATAC2](https://github.com/kaizhang/SnapATAC2) | Data frame containing the chromosome sizes for the genome build (i.e., hg38); created using the [`chrom_sizes` pipeline input](README.md). |
| `reference_gtf_file` | N/A | Path of the GTF used for TSS enrichment. |
| `library_tsse`, `frac_overlap_TSS`, `TSS_profile` | [SnapATAC2](https://github.com/kaizhang/SnapATAC2) | Library-level TSS enrichment score, fraction of fragments overlapping TSSs, and the aggregate TSS profile, computed over all barcodes. |
| `macs3_pseudobulk` | [SnapATAC2](https://github.com/kaizhang/SnapATAC2) | Peaks called by MACS3 on the pseudobulk of all barcodes. |
| `NHashID` | N/A | A string that represents the NHashID if specified in the workflow |


## Table 2. Cell metrics

| Cell Metrics | Program | Details |
|---|---|--------------------|
| `tsse` | [SnapATAC2](https://github.com/kaizhang/SnapATAC2) | Transcription start site enrichment (TSSe) score; lower scores suggest poor data quality. Learn more about TSSe in the [Definitions section](#definitions) below. |
| `n_fragment` | [SnapATAC2](https://github.com/kaizhang/SnapATAC2) | Number of unique fragments corresponding to the ATAC cell barcode. Fragments are stored in the `h5ad.obsm` property of the output h5ad file. Learn more about cell barcodes and fragments in the [Definitions section](#definitions) below. |
| `frac_dup` | [SnapATAC2](https://github.com/kaizhang/SnapATAC2) | Fraction of reads associated with the cell barcode that are duplicates. |
| `frac_mito` | [SnapATAC2](https://github.com/kaizhang/SnapATAC2) | Fraction of reads associated with the cell barcode that are mitochondrial. |
| `n_frag_overlap_peak` | [SnapATAC2](https://github.com/kaizhang/SnapATAC2) | Number of unique fragments overlapping the MACS3 pseudobulk peaks. This is the quantity used for cell calling and for the filtered h5ad ranking. |
| `is_cell` | [SnapATAC2](https://github.com/kaizhang/SnapATAC2) | Boolean; true for barcodes called as cells by the anchored OrdMag caller (counted in the `number_of_cells` library metric). |
| `is_retained` | WARP | Boolean; true for the top `atac_expected_cells` barcodes by `n_frag_overlap_peak`, i.e. the rows present in the filtered h5ad. All false when `filter_cells_to_expected` is false. |
| `gex_barcodes` | [AnnData](https://anndata.readthedocs.io/en/latest/index.html) | Gene expression barcode associated with each ATAC cell barcode. This column is only produced when ATAC is run as part of the Multiome pipeline. |


## Definitions
* Cell Barcode: Short nucleotide sequence used to label and distinguish which reads come from each unique cell, allowing for tracking of many cells simultaneously.
* Fragment: A distinct segment of a read that aligns to a specific location on the reference genome. 
* Transcription Start Site Enrichment (TSSe): A common quality control metric in ATAC-seq data, indicating increased accessibility around the transcription start sites of genes. High TSSe suggests successful capture of relevant genomic features, while low TSSe may signal data quality or processing issues.
