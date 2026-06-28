# scANVI Pipeline

## Documentation

Full scANVI pipeline documentation has moved to the [WARP documentation site](https://broadinstitute.github.io/warp/docs/Pipelines/scANVI_Pipeline/README).

## Summary

scANVI is a cloud-optimized WDL workflow that performs **cell type label transfer** using SCVI and SCANVI deep generative models. It integrates single-cell RNA-seq (GEX) and (optionally) ATAC-seq data with an annotated reference and transfers cell type labels via semi-supervised learning.

ATAC is optional. When no ATAC h5ad is supplied, the pipeline auto-detects **GEX-only mode** and trains/annotates from the reference atlas using gene expression and reference data alone, without using ATAC.

The pipeline is split into two tasks:

- **PreprocessFilter** (CPU-only) — loads and filters the input h5ad files. In multiome mode it also aligns GEX/ATAC barcodes and converts the ATAC cell-by-bin matrix to a gene activity matrix; in GEX-only mode it skips all ATAC steps.
- **MultiomeLabelTransfer** (GPU) — trains SCVI / SCANVI models, transfers labels, and writes annotated outputs. Uses `run_multi_model` (GEX + ATAC + reference) or `run_gex_only_model` (GEX + reference) depending on whether ATAC is present.

## Running the pipeline

scANVI can be run with [Cromwell](https://cromwell.readthedocs.io/en/stable/) or in [Terra](https://app.terra.bio).

Example input JSON files are in [`example_inputs/`](./example_inputs).

Required inputs:

- `scANVI.input_id` — unique identifier prepended to all output filenames.
- Either `scANVI.input_bucket` (a GCS path containing `gex.h5ad`, `ref.h5ad`, and optionally `atac.h5ad`) **or** the direct file inputs `scANVI.gex_h5ad` and `scANVI.ref_h5ad`, plus optionally `scANVI.atac_h5ad`.

The ATAC input (`scANVI.atac_h5ad`, or `atac.h5ad` in the bucket) is **optional**. If it is omitted/absent, the pipeline runs in GEX-only mode and the `atac_annotated_h5ad` output is not produced. See [`example_inputs/scANVI.gex_only.json`](./example_inputs/scANVI.gex_only.json) for a GEX-only example.

The reference h5ad must contain cell type annotations in `obs['final_annotation']`.

## Versioning

See [scANVI.changelog.md](scANVI.changelog.md) for the full release history.
