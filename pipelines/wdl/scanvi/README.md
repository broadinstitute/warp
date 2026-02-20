# scANVI (ScviScanvi) Pipeline

## Overview

The ScviScanvi pipeline is a cloud-optimized WDL workflow for performing **cell type label transfer on Multiome data** using [scVI](https://docs.scvi-tools.org/en/1.4.1/user_guide/models/scvi.html) and [scANVI](https://docs.scvi-tools.org/en/1.4.1/user_guide/models/scanvi.html) (single-cell ANnotation using Variational Inference) deep generative models. It integrates single-cell RNA (GEX) and ATAC data with an annotated reference dataset to transfer cell type labels via semi-supervised learning.

For detailed information about the scANVI model, see the [scANVI documentation](https://docs.scvi-tools.org/en/1.4.1/user_guide/models/scanvi.html). For details on the Docker image and underlying scripts, see the [warp-tools scvi-scanvi README](https://github.com/broadinstitute/warp-tools/tree/develop/3rd-party-tools/scvi-scanvi).

## Inputs

| Input | Type | Description | Default |
|---|---|---|---|
| `gex_h5ad` | File | Gene expression AnnData h5ad file from Multiome/Optimus pipeline output | Required |
| `atac_h5ad` | File | ATAC cell-by-bin AnnData h5ad file from Multiome/PeakCalling pipeline output | Required |
| `ref_h5ad` | File | Annotated reference AnnData h5ad file with cell type labels in `obs['final_annotation']` | Required |
| `cloud_provider` | String | Cloud platform: `"gcp"` or `"azure"` | Required |
| `disk_size` | Int | Disk size in GB | 500 |
| `mem_size` | Int | Memory size in GB | 64 |
| `nthreads` | Int | Number of CPU threads | 8 |
| `gpu_type` | String | GPU type for accelerated model training | `"nvidia-tesla-t4"` |
| `gpu_count` | Int | Number of GPUs | 1 |

## Outputs

| Output | Type | Description |
|---|---|---|
| `scanvi_predictions_h5ad` | File | SCANVI cell type predictions as an h5ad file |
| `gex_annotated_h5ad` | File | Gene expression AnnData annotated with transferred cell type labels |
| `atac_annotated_h5ad` | File | ATAC AnnData annotated with transferred cell type labels |
| `pipeline_version_out` | String | Pipeline version string |

## How It Works

The pipeline runs a single task (`MultiomeLabelTransfer`) that performs the following steps:

1. **Preprocessing** — Filters GEX data, reindexes ATAC barcodes to align with GEX, retains shared barcodes, and aligns gene features across GEX, ATAC, and the reference dataset.
2. **SCVI training** — Trains an scVI model on the combined data to learn an unsupervised latent representation.
3. **SCANVI training** — Extends the scVI model with semi-supervised learning, using cell type annotations from the reference to transfer labels to unlabelled query cells.
4. **Label transfer** — Applies the trained SCANVI model to predict cell types and produces annotated GEX and ATAC outputs.

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
