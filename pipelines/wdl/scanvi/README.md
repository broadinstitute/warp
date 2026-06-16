# scANVI Pipeline

## Documentation

Full scANVI pipeline documentation has moved to the [WARP documentation site](https://broadinstitute.github.io/warp/docs/Pipelines/scANVI_Pipeline/README).

## Summary

scANVI is a cloud-optimized WDL workflow that performs **cell type label transfer on Multiome data** using SCVI and SCANVI deep generative models. It integrates single-cell RNA-seq (GEX) and ATAC-seq data with an annotated reference and transfers cell type labels via semi-supervised learning.

The pipeline is split into two tasks:

- **PreprocessFilter** (CPU-only) — loads, filters, and aligns the three input h5ad files; converts the ATAC cell-by-bin matrix to a gene activity matrix.
- **MultiomeLabelTransfer** (GPU) — trains SCVI / SCANVI models, transfers labels, and writes annotated outputs.

## Running the pipeline

scANVI can be run with [Cromwell](https://cromwell.readthedocs.io/en/stable/) or in [Terra](https://app.terra.bio).

Example input JSON files are in [`example_inputs/`](./example_inputs).

Required inputs:

- `scANVI.input_id` — unique identifier prepended to all output filenames.
- Either `scANVI.input_bucket` (a GCS path containing `gex.h5ad`, `atac.h5ad`, `ref.h5ad`) **or** the three direct file inputs `scANVI.gex_h5ad`, `scANVI.atac_h5ad`, `scANVI.ref_h5ad`.

The reference h5ad must contain cell type annotations in `obs['final_annotation']`.

## Versioning

See [scANVI.changelog.md](scANVI.changelog.md) for the full release history.
