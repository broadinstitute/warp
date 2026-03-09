# scANVI Pipeline

## Overview

The scANVI pipeline is a cloud-optimized WDL workflow that performs **cell type label transfer on Multiome data** using [scVI](https://docs.scvi-tools.org/en/1.4.1/user_guide/models/scvi.html) and [scANVI](https://docs.scvi-tools.org/en/1.4.1/user_guide/models/scanvi.html) (single-cell ANnotation using Variational Inference) deep generative models. It integrates single-cell RNA-seq (GEX) and ATAC-seq data with an annotated reference dataset to transfer cell type labels via semi-supervised learning.

The pipeline is split into two tasks — **PreprocessFilter** (CPU-only) and **MultiomeLabelTransfer** (GPU) — so that expensive GPU time is reserved exclusively for model training and inference. All data loading, quality filtering, barcode alignment, and gene-activity-matrix conversion happen on a CPU node first; the GPU node receives ready-to-train h5ad files and never re-runs preprocessing.

## Inputs

| Input | Type | Description | Default |
|---|---|---|---|
| `input_bucket` | String? | GCS bucket path containing input h5ad files | — |
| `gex_h5ad` | File? | Gene expression h5ad from Multiome/Optimus | — |
| `atac_h5ad` | File? | ATAC cell-by-bin h5ad from Multiome/PeakCalling | — |
| `ref_h5ad` | File? | Annotated reference h5ad with `obs['final_annotation']` | — |
| `gex_filename` | String | GEX filename in the bucket | `"gex.h5ad"` |
| `atac_filename` | String | ATAC filename in the bucket | `"atac.h5ad"` |
| `ref_filename` | String | Reference filename in the bucket | `"ref.h5ad"` |
| `disk_size` | Int | Disk size (GB) | 500 |
| `mem_size` | Int | Memory (GB) | 120 |
| `nthreads` | Int | CPU threads | 32 |

**Input modes:** Provide either `input_bucket` (to download all three h5ad files) or individual `File` inputs (`gex_h5ad`, `atac_h5ad`, `ref_h5ad`). Direct file inputs take precedence.

## Outputs

| Output | Type | Description |
|---|---|---|
| `scanvi_predictions_h5ad` | File | Combined AnnData with SCANVI cell type predictions, UMAP, and metadata |
| `gex_annotated_h5ad` | File | Preprocessed GEX AnnData annotated with transferred cell type labels |
| `atac_annotated_h5ad` | File | ATAC gene-activity AnnData annotated with transferred cell type labels |
| `pipeline_version_out` | String | Pipeline version string |

## How It Works

### Task 1 — `PreprocessFilter` (CPU-only)

Loads and preprocesses the three input h5ad files on a CPU node. No GPU is allocated for this task. Steps:

1. **Load datasets** — Reads GEX (scanpy), ATAC cell-by-bin (snapatac2), and reference (scanpy).
2. **Patch missing columns** — Adds `star_IsCell = True` to GEX and `gex_barcodes` (from index) to ATAC if absent, ensuring compatibility across upstream pipelines.
3. **Filter GEX** — Retains STARsolo cell calls (`star_IsCell == True`), then removes genes and cells with < 3 counts.
4. **Prepare GEX** — Sets `batch` label; copies counts into a `counts` layer.
5. **Reindex ATAC** — Sets ATAC obs index to `gex_barcodes` so barcodes align with GEX.
6. **Shared barcode filtering** — Intersects GEX and ATAC barcodes; subsets both to matched cells.
7. **Batch labels** — GEX → `pd-multiome_sci_gex`, ATAC → `pd-multiome_sci_atac`.
8. **Placeholder annotations** — Adds `final_annotation = "Unknown"` to query datasets.
9. **Gene activity matrix** — Converts the ATAC cell-by-bin matrix into a gene activity matrix via `snapatac2.pp.make_gene_matrix` (hg38 GENCODE annotation).
10. **Modality tags** — GEX → `rna_unannotated`, ATAC activity → `atac_unannotated`, reference → `rna_annotated`.
11. **Write outputs** — `preprocessed_gex.h5ad`, `preprocessed_atac_activity.h5ad`, `preprocessed_ref.h5ad`.

### Task 2 — `MultiomeLabelTransfer` (GPU)

Loads the three preprocessed h5ad files and performs **only** model training, label transfer, and output finalization. It imports individual functions (`run_multi_model`, `transfer_labels`, `finalize_output`) from the container's `multiome_label_transfer.py` module — the script's `main()` function is **never called**, so no preprocessing is repeated.

1. **Load preprocessed data** — Reads the three h5ad files produced by PreprocessFilter. No filtering, reindexing, or conversion is performed.
2. **Train SCVI** — `run_multi_model()` concatenates the three AnnData objects, filters to genes in ≥ 5 cells, selects 5000 highly variable genes (Seurat v3, batch-aware), then trains an **SCVI** model (unsupervised VAE: 2 layers, 30 latent dimensions, negative-binomial likelihood, gene-batch dispersion, up to 500 epochs with early stopping).
3. **Train SCANVI** — The same function initializes **SCANVI** from the trained SCVI model and performs semi-supervised training using the reference cell type labels (`final_annotation`), propagating annotations to unlabeled GEX and ATAC cells (up to 500 epochs, 100 samples per label).
4. **Predict labels** — `transfer_labels()` uses the trained SCANVI model to predict cell types (`C_scANVI`) for every cell, extracts the latent representation (`X_scANVI`), computes a neighborhood graph and UMAP embedding.
5. **Propagate labels** — Copies the predicted `C_scANVI` labels from the concatenated object back into the original GEX and ATAC AnnData objects using the barcode-suffix index created by `ad.concat`.
6. **Write annotated matrices** — Saves `gex_annotated_matrix.h5ad` and `atac_annotated_matrix.h5ad`.
7. **Finalize predictions** — `finalize_output()` adds placeholder metadata (biosample, donor, species, disease, organ, library prep, sex), renames `final_annotation` → `celltype`, and copies counts into the `.raw` layer for SCP ingest. Writes `SCANVI_predictions.h5ad`.

### Workflow Diagram

```
                ┌──────────────────────┐
                │   Input h5ad files   │
                │  (GEX, ATAC, Ref)    │
                └──────────┬───────────┘
                           │
                           ▼
              ┌────────────────────────┐
              │   PreprocessFilter     │  CPU-only
              │  (load, filter, align, │
              │   gene activity matrix)│
              └──────────┬─────────────┘
                         │
          ┌──────────────┼──────────────┐
          ▼              ▼              ▼
   preprocessed    preprocessed    preprocessed
     _gex.h5ad    _atac_activity     _ref.h5ad
                      .h5ad
          │              │              │
          └──────────────┼──────────────┘
                         │
                         ▼
            ┌────────────────────────────┐
            │  MultiomeLabelTransfer     │  GPU
            │  (import run_multi_model,  │
            │   transfer_labels,         │
            │   finalize_output —        │
            │   main() is NOT called)    │
            └──────────┬─────────────────┘
                       │
          ┌────────────┼────────────┐
          ▼            ▼            ▼
   SCANVI_predictions  gex_annotated  atac_annotated
        .h5ad          _matrix.h5ad   _matrix.h5ad
```

## Design Rationale

The container script `multiome_label_transfer.py` bundles preprocessing **and** model training inside a single `main()` function. If the GPU task called `main()`, all preprocessing would run a second time on the GPU node — wasting expensive GPU hours on CPU-bound work that has already been completed.

Instead, the pipeline:
- Runs all preprocessing in **PreprocessFilter** on a CPU-only VM (no GPU cost).
- In **MultiomeLabelTransfer**, imports only the three model-training / label-transfer / finalization functions from the script, bypassing `main()` entirely.

This ensures zero duplication: every preprocessing step executes exactly once (on CPU), and the GPU node spends 100 % of its time on model training and inference.

## Requirements

- **GPU required** for Task 2 — SCVI/SCANVI training benefits significantly from GPU acceleration (NVIDIA Tesla T4 or equivalent).
- Input h5ad files should follow [Multiome](https://broadinstitute.github.io/warp/docs/Pipelines/Multiome_Pipeline/README) and [PeakCalling](https://broadinstitute.github.io/warp/docs/Pipelines/PeakCalling/README) WARP pipeline conventions.
- The reference h5ad must contain cell type annotations in `obs['final_annotation']`.

## Docker

Both tasks use the same `scvi-scanvi` Docker image from [warp-tools](https://github.com/broadinstitute/warp-tools/tree/develop/3rd-party-tools/scvi-scanvi). Only Task 2 gets GPUs attached by the execution engine (hardcoded as `nvidia-tesla-t4` x 2 in the runtime block); the container itself does not configure CUDA — the execution engine handles GPU/driver setup.

```
us.gcr.io/broad-gotc-prod/scvi-scanvi:1.0.0-1.2-1756234975
```

Key libraries: scvi-tools 1.2, snapatac2 2.7, scanpy, anndata.

## Versioning and Changelog

See [scANVI.changelog.md](scANVI.changelog.md) for the full release history.
