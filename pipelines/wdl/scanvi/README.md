# scANVI Pipeline

## Documentation

Full scANVI pipeline documentation has moved to the [WARP documentation site](https://broadinstitute.github.io/warp/docs/Pipelines/scANVI_Pipeline/README).

## Summary

scANVI is a cloud-optimized WDL workflow that performs **cell type label transfer** using SCVI and SCANVI deep generative models. It integrates single-cell RNA-seq (GEX) and (optionally) ATAC-seq data with an annotated reference and transfers cell type labels via semi-supervised learning.

ATAC is optional. When no ATAC h5ads are supplied, the pipeline runs in **GEX-only mode** and trains/annotates from the reference atlas using gene expression and reference data alone, without using ATAC.

One run labels one or more samples, each independently and with its own outputs, on one machine, one sample after another. Each sample is labelled at one or more label levels (for example `class` and `subclass`), one SCANVI model per level, and each input file `<name>.h5ad` gives one `<name>_annotated.h5ad`, written ready for Single Cell Portal (counts in `.raw`, normalized `.X`, and per level the labels, a 2-D UMAP and the scVI latent). Every output records its sample's id, its input file name, its model's name, which is the reference the model was trained on, and the levels in `uns`. No sample metadata is added.

The pipeline has two tasks:

- **PreprocessReference** (CPU-only, once per run, training only) — adapts the annotated reference (AIT or PBMC-style) into the form the trainer expects.
- **MultiomeLabelTransfer** (GPU, or CPU as MultiomeLabelTransferCpu) — checks every input, then preprocesses and labels each sample in turn: trains SCVI and a SCANVI per level on the reference plus the sample, or uses a supplied model archive, mapping the sample onto it when the model was trained on the reference alone.

## Running the pipeline

scANVI can be run with [Cromwell](https://cromwell.readthedocs.io/en/stable/) or in [Terra](https://app.terra.bio).

Example input JSON files are in [`example_inputs/`](./example_inputs).

Required inputs:

- `scANVI.input_ids` and `scANVI.gex_h5ads` — arrays paired by position, one entry per sample; `input_ids[i]` is recorded in every output of `gex_h5ads[i]`. Outputs are named after the input files, which must have distinct names. A single sample is a one-element array (in JSON, `["gs://..."]`).
- `scANVI.scanvi_model` (a model archive, used for every sample) or `scANVI.ref_h5ad` (the annotated reference, to train on).

The ATAC input `scANVI.atac_h5ads` is **optional**: leave it empty (the default) for GEX-only mode, or give one ATAC h5ad per `gex_h5ads` entry, in the same order, for multiome mode. In GEX-only mode `atac_annotated_h5ads` is empty. See [`example_inputs/scANVI.gex_only.json`](./example_inputs/scANVI.gex_only.json) for a GEX-only example.

A supplied `scANVI.scanvi_model` is used for every sample and the reference is not needed; without one, each sample trains its own model, so training cost grows with the number of samples.

When training, the reference h5ad must provide cell-type labels and a batch. **PBMC-style** references have counts in `.X`, labels in `obs['final_annotation']`, and a batch in `obs['batch']`. **AIT-schema** (Allen Institute Taxonomy) references are auto-detected and adapted: counts are pulled from `.raw`, and the label/batch default to `subclass`/`donor_id` (override with `scANVI.ref_label_columns`, one or more label levels, and `scANVI.ref_batch_column`). For a mouse AIT reference in multiome mode, set `scANVI.genome` to `mm10`/`mm39` (default `hg38`).

### Per-cell label confidence

Every annotated h5ad has a `<level>_max_probability` obs column per level — the per-cell maximum SCANVI posterior probability, i.e. the confidence of the assigned label, useful for thresholding or flagging low-confidence calls downstream. Set `scANVI.output_max_probability` to `false` to leave it out.

## Versioning

See [scANVI.changelog.md](scANVI.changelog.md) for the full release history.

## Citing

If you use the scANVI Pipeline in your research, please identify the pipeline in your methods section using the [scANVI SciCrunch resource identifier](https://rrid.site/data/record/nlx_144509-1/SCR_028705/resolver?q=scANVI&l=scANVI&i=rrid:scr_028705).

* Ex: *scANVI Pipeline (RRID:SCR_028705)*

See [Citing the scANVI Pipeline](https://broadinstitute.github.io/warp/docs/Pipelines/scANVI_Pipeline/README#citing-the-scanvi-pipeline) for the full WARP citation.
