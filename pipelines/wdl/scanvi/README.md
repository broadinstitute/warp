# scANVI Pipeline

## Documentation

Full scANVI pipeline documentation has moved to the [WARP documentation site](https://broadinstitute.github.io/warp/docs/Pipelines/scANVI_Pipeline/README).

## Summary

scANVI is a cloud-optimized WDL workflow that performs **cell type label transfer** using SCVI and SCANVI deep generative models. It integrates single-cell RNA-seq (GEX) and (optionally) ATAC-seq data with an annotated reference and transfers cell type labels via semi-supervised learning.

ATAC is optional. When no ATAC h5ads are supplied, the pipeline runs in **GEX-only mode** and trains/annotates from the reference atlas using gene expression and reference data alone, without using ATAC.

One run labels one or more query files against a single reference. The workflow scatters over the queries: each is labelled independently and gets its own outputs.

The pipeline is split into three tasks:

- **PreprocessReference** (CPU-only, once per run) — adapts the annotated reference (AIT or PBMC-style) into the form the trainer expects.
- **PreprocessFilter** (CPU-only, once per query) — loads and filters one query's GEX h5ad. In multiome mode it also aligns GEX/ATAC barcodes and converts the ATAC cell-by-bin matrix to a gene activity matrix; in GEX-only mode it skips all ATAC steps.
- **MultiomeLabelTransfer** (GPU, once per query) — trains SCVI / SCANVI models (or loads a supplied one), transfers labels, and writes annotated outputs. Uses `run_multi_model` (GEX + ATAC + reference) or `run_gex_only_model` (GEX + reference) depending on whether ATAC is present.

## Running the pipeline

scANVI can be run with [Cromwell](https://cromwell.readthedocs.io/en/stable/) or in [Terra](https://app.terra.bio).

Example input JSON files are in [`example_inputs/`](./example_inputs).

Required inputs:

- `scANVI.input_ids` and `scANVI.gex_h5ads` — arrays paired by position, one entry per query; `input_ids[i]` prefixes every output of `gex_h5ads[i]`. A single query is a one-element array (in JSON, `["gs://..."]`).
- `scANVI.ref_h5ad` — the annotated reference, shared by every query.

The ATAC input `scANVI.atac_h5ads` is **optional**: leave it empty (the default) for GEX-only mode, or give one ATAC h5ad per `gex_h5ads` entry, in the same order, for multiome mode. In GEX-only mode `atac_annotated_h5ads` is empty. See [`example_inputs/scANVI.gex_only.json`](./example_inputs/scANVI.gex_only.json) for a GEX-only example.

A supplied `scANVI.scanvi_model` is used for every query; without one, each query trains its own model, so training cost grows with the number of queries.

The reference h5ad must provide cell-type labels and a batch. **PBMC-style** references have counts in `.X`, labels in `obs['final_annotation']`, and a batch in `obs['batch']`. **AIT-schema** (Allen Institute Taxonomy) references are auto-detected and adapted: counts are pulled from `.raw`, and the label/batch default to `subclass`/`donor_id` (override with `scANVI.ref_label_column` / `scANVI.ref_batch_column`). For a mouse AIT reference in multiome mode, set `scANVI.genome` to `mm10`/`mm39` (default `hg38`).

### Optional: per-cell label confidence

Set `scANVI.output_max_probability` to `true` to add a `max_probability` obs column to every output h5ad — the per-cell maximum SCANVI posterior probability, i.e. the confidence of the assigned label. Default `false` leaves outputs unchanged. Useful for thresholding or flagging low-confidence calls downstream.

## Versioning

See [scANVI.changelog.md](scANVI.changelog.md) for the full release history.

## Citing

If you use the scANVI Pipeline in your research, please identify the pipeline in your methods section using the [scANVI SciCrunch resource identifier](https://rrid.site/data/record/nlx_144509-1/SCR_028705/resolver?q=scANVI&l=scANVI&i=rrid:scr_028705).

* Ex: *scANVI Pipeline (RRID:SCR_028705)*

See [Citing the scANVI Pipeline](https://broadinstitute.github.io/warp/docs/Pipelines/scANVI_Pipeline/README#citing-the-scanvi-pipeline) for the full WARP citation.
