# scANVI (ScviScanvi) Pipeline

## Overview

The ScviScanvi pipeline is a cloud-optimized WDL workflow for performing **cell type label transfer on Multiome data** using [scVI](https://docs.scvi-tools.org/en/1.4.1/user_guide/models/scvi.html) and [scANVI](https://docs.scvi-tools.org/en/1.4.1/user_guide/models/scanvi.html) (single-cell ANnotation using Variational Inference) deep generative models. It integrates single-cell RNA (GEX) and ATAC data with an annotated reference dataset to transfer cell type labels via semi-supervised learning.

For detailed information about the scANVI model, see the [scANVI documentation](https://docs.scvi-tools.org/en/1.4.1/user_guide/models/scanvi.html). For details on the Docker image and underlying scripts, see the [warp-tools scvi-scanvi README](https://github.com/broadinstitute/warp-tools/tree/develop/3rd-party-tools/scvi-scanvi).

## Inputs

| Input | Type | Description | Default |
|---|---|---|---|
| `input_bucket` | String? | GCS bucket path containing input h5ad files (e.g., `gs://bucket/path/to/inputs`) | Optional |
| `gex_h5ad` | File? | Gene expression AnnData h5ad file from Multiome/Optimus pipeline output | Optional |
| `atac_h5ad` | File? | ATAC cell-by-bin AnnData h5ad file from Multiome/PeakCalling pipeline output | Optional |
| `ref_h5ad` | File? | Annotated reference AnnData h5ad file with cell type labels in `obs['final_annotation']` | Optional |
| `gex_filename` | String | Expected GEX h5ad filename in the input bucket | `"gex.h5ad"` |
| `atac_filename` | String | Expected ATAC h5ad filename in the input bucket | `"atac.h5ad"` |
| `ref_filename` | String | Expected reference h5ad filename in the input bucket | `"ref.h5ad"` |
| `cloud_provider` | String | Cloud platform | `"gcp"` |
| `disk_size` | Int | Disk size in GB | 500 |
| `mem_size` | Int | Memory size in GB | 120 |
| `nthreads` | Int | Number of CPU threads | 32 |
| `gpu_type` | String | GPU type for accelerated model training | `"nvidia-tesla-t4"` |
| `gpu_count` | Int | Number of GPUs | 2 |
| `nvidiaDriverVersion` | String | NVIDIA driver version for GPU support | `"535.104.05"` |

**Input modes:** Provide either `input_bucket` (to download all three h5ad files from a GCS bucket) or individual `File` inputs (`gex_h5ad`, `atac_h5ad`, `ref_h5ad`). Direct file inputs take precedence over bucket-downloaded files.

## Outputs

| Output | Type | Description |
|---|---|---|
| `scanvi_predictions_h5ad` | File | SCANVI cell type predictions as an h5ad file |
| `gex_annotated_h5ad` | File | Gene expression AnnData annotated with transferred cell type labels |
| `atac_annotated_h5ad` | File | ATAC AnnData annotated with transferred cell type labels |
| `pipeline_version_out` | String | Pipeline version string |

## How It Works

The pipeline is split into two tasks to separate CPU-only preprocessing from GPU-accelerated model training:

### Task 1: `PreprocessFilter` (CPU-only)

Handles all h5ad preprocessing and filtering before model training:

1. **Load** — Reads the GEX (scanpy), ATAC cell-by-bin (snapatac2), and reference (scanpy) h5ad files.
2. **Patch missing columns** — Adds `star_IsCell` (all `True`) to GEX and `gex_barcodes` (copy of obs index) to ATAC if the columns are absent, ensuring compatibility with inputs from different upstream pipelines.
3. **Filter GEX cells** — Retains only barcodes where `star_IsCell == True` (STARsolo cell calls), then removes genes and cells with fewer than 3 total counts (`scanpy.pp.filter_genes`, `scanpy.pp.filter_cells`).
4. **Prepare GEX** — Adds a `batch` column and copies the count matrix into a `counts` layer.
5. **Reindex ATAC barcodes** — Sets the ATAC obs index to the `gex_barcodes` column so barcodes align with the GEX dataset.
6. **Shared barcode filtering** — Intersects GEX and ATAC barcodes and subsets both to only matched cells.
7. **Assign batch labels & modality tags** — Labels GEX as `pd-multiome_sci_gex` / `rna_unannotated`, ATAC as `pd-multiome_sci_atac` / `atac_unannotated`, reference as `rna_annotated`.
8. **Convert ATAC to gene activity** — Transforms the ATAC cell-by-bin matrix into a gene activity matrix using `snapatac2.pp.make_gene_matrix` with the hg38 GENCODE annotation.
9. **Set placeholder annotations** — Adds `final_annotation = "Unknown"` to unannotated query datasets.
10. **Write outputs** — Produces `preprocessed_gex.h5ad`, `preprocessed_atac_activity.h5ad`, and `preprocessed_ref.h5ad`.

### Task 2: `MultiomeLabelTransfer` (GPU)

Receives the three preprocessed h5ad files and runs model training and label transfer via `multiome_label_transfer.py`:

1. **SCVI training** — Trains an scVI model on the combined data to learn an unsupervised latent representation.
2. **SCANVI training** — Extends the scVI model with semi-supervised learning, using cell type annotations from the reference to transfer labels to unlabelled query cells.
3. **Label transfer** — Applies the trained SCANVI model to predict cell types and produces annotated GEX and ATAC outputs.

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
              │  (patch, filter, align,│
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
            ┌────────────────────────┐
            │ MultiomeLabelTransfer  │  GPU
            │  (SCVI → SCANVI →     │
            │   label transfer)     │
            └──────────┬─────────────┘
                       │
          ┌────────────┼────────────┐
          ▼            ▼            ▼
   SCANVI_predictions  gex_annotated  atac_annotated
        .h5ad          _matrix.h5ad   _matrix.h5ad
```

For a full description of the model architecture and training procedure, see the [scANVI user guide](https://docs.scvi-tools.org/en/1.4.1/user_guide/models/scanvi.html).

## Requirements

- **GPU strongly recommended** — scANVI model training benefits significantly from GPU acceleration (NVIDIA Tesla T4 or equivalent).
- Input h5ad files are expected to follow the conventions of the [Multiome](https://broadinstitute.github.io/warp/docs/Pipelines/Multiome_Pipeline/README) and [PeakCalling](https://broadinstitute.github.io/warp/docs/Pipelines/PeakCalling/README) WARP pipelines.
- The reference h5ad must contain cell type annotations in `obs['final_annotation']`.

## Docker

The pipeline uses the `scvi-scanvi` Docker image maintained in [warp-tools](https://github.com/broadinstitute/warp-tools/tree/develop/3rd-party-tools/scvi-scanvi):

```
us.gcr.io/broad-gotc-prod/scvi-scanvi:1.0.0-1.2-1756234975
```

## Example Inputs

An example input JSON is available at [`example_inputs/ScviScanvi.json`](example_inputs/ScviScanvi.json).

## Versioning and Changelog

See [scANVI.changelog.md](scANVI.changelog.md) for the full release history.
