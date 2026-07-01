---
sidebar_position: 1
slug: /Pipelines/scANVI_Pipeline/README
---

# scANVI Overview

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [scANVI_v1.2.0](https://github.com/broadinstitute/warp/releases?q=scANVI&expanded=true) | June, 2026 | WARP Pipelines | Please [file an issue in WARP](https://github.com/broadinstitute/warp/issues) |

## Introduction to the scANVI workflow

The scANVI pipeline is a cloud-optimized WDL workflow that performs **cell type label transfer** using [scVI](https://docs.scvi-tools.org/en/1.4.1/user_guide/models/scvi.html) and [scANVI](https://docs.scvi-tools.org/en/1.4.1/user_guide/models/scanvi.html) (single-cell ANnotation using Variational Inference) deep generative models. It integrates single-cell RNA-seq (GEX) and, optionally, ATAC-seq data with an annotated reference dataset to transfer cell type labels via semi-supervised learning.

ATAC is optional. When no ATAC h5ad is supplied, the pipeline auto-detects **GEX-only mode** and trains/annotates from the reference atlas using gene expression and reference data alone, without using ATAC. When ATAC is supplied it runs in multiome mode (GEX + ATAC + reference).

The pipeline is split into two tasks — **PreprocessFilter** (CPU-only) and **MultiomeLabelTransfer** (GPU) — so that expensive GPU time is reserved exclusively for model training and inference. All data loading, quality filtering, and (in multiome mode) barcode alignment and gene-activity-matrix conversion happen on a CPU node first; the GPU node receives ready-to-train h5ad files and never re-runs preprocessing.

:::tip Want to use scANVI for your publication?
The pipeline is designed to consume the outputs of the [Multiome](../Multiome_Pipeline/README.md) and [PeakCalling](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/peak_calling) WARP pipelines. Cite the pipeline using the WARP citation in the [Citing](#citing-the-scanvi-pipeline) section below.
:::

## Quickstart table

The following table provides a quick glance at the scANVI pipeline features:

| Pipeline features | Description | Source |
| --- | --- | --- |
| Assay type | 10x single-cell / single-nucleus Multiome (GEX + ATAC), or GEX-only | [10x Genomics](https://www.10xgenomics.com) |
| Overall workflow | CPU preprocessing + GPU SCVI/SCANVI label transfer | Code available on [GitHub](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.wdl) |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Models | SCVI (unsupervised VAE) + SCANVI (semi-supervised classifier) | [scvi-tools 1.2](https://docs.scvi-tools.org/) |
| ATAC gene-activity conversion | Cell-by-bin matrix → gene activity matrix (hg38 GENCODE) | [snapatac2 2.7](https://kzhang.org/SnapATAC2/) |
| Data input format | AnnData h5ad files: GEX and annotated reference (required), ATAC cell-by-bin (optional) | [AnnData](https://anndata.readthedocs.io/) |
| Data output format | Annotated h5ad files with predicted cell types and UMAP | [AnnData](https://anndata.readthedocs.io/) |

## Set-up

### scANVI installation

To download the latest scANVI release, see the release tags prefixed with "scANVI" on the WARP [releases page](https://github.com/broadinstitute/warp/releases). All scANVI pipeline releases are documented in the [scANVI changelog](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.changelog.md).

scANVI can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant, flexible workflow management system that supports multiple computing platforms. The workflow can also be run in [Terra](https://app.terra.bio).

### Inputs

scANVI accepts inputs in two modes — **direct file inputs** or **bucket mode**. Direct file inputs take precedence when both are provided.

Example input JSON files are available in the [`example_inputs`](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/scanvi/example_inputs) folder.

#### Workflow inputs

| Input | Type | Description | Default |
| --- | --- | --- | --- |
| `input_id` | String | Unique identifier prepended to all output filenames. | — (required) |
| `input_bucket` | String? | GCS bucket path containing input h5ad files (e.g., `gs://bucket/path/to/inputs`). Used when direct file inputs are not provided. | — |
| `gex_h5ad` | File? | Gene expression AnnData h5ad file from Multiome / Optimus. | — |
| `atac_h5ad` | File? | **Optional** ATAC cell-by-bin AnnData h5ad file from Multiome / PeakCalling. If omitted, the pipeline runs in GEX-only mode. | — |
| `ref_h5ad` | File? | Annotated reference AnnData h5ad file with cell type labels in `obs['final_annotation']`. | — |
| `gex_filename` | String | Expected GEX h5ad filename in the input bucket. | `"gex.h5ad"` |
| `atac_filename` | String | Expected ATAC h5ad filename in the input bucket. Optional: if absent from the bucket, the pipeline runs in GEX-only mode. | `"atac.h5ad"` |
| `ref_filename` | String | Expected reference h5ad filename in the input bucket. | `"ref.h5ad"` |
| `max_epochs` | Int? | Optional cap on SCVI/SCANVI training epochs, applied in both multiome and GEX-only modes. When unset, the container default (500) is used. | — |
| `ref_label_column` | String? | Reference `obs` column to use as the cell-type label. When unset, defaults to `subclass` for AIT references and `final_annotation` otherwise. | — |
| `ref_batch_column` | String? | Reference `obs` column to use as the batch. When unset, defaults to `donor_id` for AIT references and `batch` otherwise. | — |
| `genome` | String | Genome for the ATAC cell-by-bin → gene-activity conversion (multiome only): `hg38` (default), `mm10`, or `mm39`. | `"hg38"` |

:::note Input mode precedence
If `gex_h5ad` and `ref_h5ad` are supplied, they are used directly and `input_bucket` is ignored. Otherwise, the filenames are downloaded from `input_bucket` via `gsutil`. The pipeline fails fast if a required (GEX or reference) input is missing or empty.
:::

:::note GEX-only mode
The ATAC input is optional. In direct-file mode, omit `atac_h5ad`; in bucket mode, simply do not place an `atac.h5ad` in the bucket (its presence is auto-detected via `gsutil stat`). In GEX-only mode, SCVI/SCANVI are trained on GEX + reference via `run_gex_only_model`, labels are transferred to GEX only, and the `atac_annotated_h5ad` output is not produced. See [`example_inputs/scANVI.gex_only.json`](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/scanvi/example_inputs).
:::

#### Reference requirements

The reference h5ad must provide cell-type labels and a batch. Two reference forms are supported:

- **PBMC-style** (default): counts in `.X`, labels in `obs['final_annotation']`, batch in `obs['batch']`.
- **AIT-schema** (Allen Institute Taxonomy): auto-detected via `uns['schema_version']` + `uns['hierarchy']`. These typically have **no `.X`** (counts are materialized from `.raw`, keeping the gene symbols from `.var`), a cell-type **hierarchy** (e.g. `class`/`subclass`/`cluster_id`), and a batch column such as `donor_id`. By default the label is taken from `subclass` and the batch from `donor_id`; override with `ref_label_column` / `ref_batch_column`. For a mouse AIT reference used in multiome mode, set `genome` to `mm10`/`mm39` so the ATAC gene-activity conversion uses the correct genome.

The query datasets (GEX, and ATAC when present) do not need pre-existing annotations — placeholder `Unknown` labels are added automatically before training.

## scANVI tasks and tools

The [scANVI workflow](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.wdl) defines two tasks inline. Both use the same Docker image; only the second task is allocated GPUs.

| Task name | Tool | Software | Description |
| --- | --- | --- | --- |
| [PreprocessFilter](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.wdl) | scanpy, snapatac2 | Python | Loads the input h5ad files, patches missing columns, and filters GEX to STARsolo cell calls. In multiome mode it also intersects barcodes between GEX and ATAC and converts the ATAC cell-by-bin matrix into a gene activity matrix; in GEX-only mode these ATAC steps are skipped. |
| [MultiomeLabelTransfer](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.wdl) | scvi-tools | Python (GPU) | Trains an SCVI model on the preprocessed AnnData objects (GEX + ATAC + reference, or GEX + reference in GEX-only mode), then trains an SCANVI classifier from the SCVI model to transfer reference cell-type labels onto the unlabeled query cells. Computes UMAP and writes annotated outputs. |

Overall, the scANVI workflow:

1. [Preprocesses and filters the three input h5ad files (CPU).](#1-preprocessfilter-cpu-only)
2. [Trains SCVI / SCANVI models and transfers labels (GPU).](#2-multiomelabeltransfer-gpu)

#### 1. PreprocessFilter (CPU-only)

Loads and preprocesses the input h5ad files on a CPU-only node. No GPU is allocated for this task. Steps 5, 6, and 9 (the ATAC-specific steps) run **only in multiome mode**; in GEX-only mode they are skipped and no ATAC activity file is written. Steps:

1. **Load datasets** — Reads GEX (scanpy) and reference (scanpy), plus ATAC cell-by-bin (snapatac2) when present.
2. **Patch missing columns** — Adds `star_IsCell = True` to GEX and `gex_barcodes` (from index) to ATAC if absent, ensuring compatibility across upstream pipelines.
3. **Filter GEX** — Retains STARsolo cell calls (`star_IsCell == True`), then removes genes and cells with fewer than 3 counts.
4. **Prepare GEX** — Sets `batch` label; copies counts into a `counts` layer.
5. **Reindex ATAC** *(multiome only)* — Sets ATAC obs index to `gex_barcodes` so barcodes align with GEX.
6. **Shared barcode filtering** *(multiome only)* — Intersects GEX and ATAC barcodes; subsets both to matched cells. In GEX-only mode all filtered GEX cells are retained.
7. **Batch labels** — GEX → `pd-multiome_sci_gex` (and, in multiome mode, ATAC → `pd-multiome_sci_atac`).
8. **Placeholder annotations** — Adds `final_annotation = "Unknown"` to query datasets.
9. **Gene activity matrix** *(multiome only)* — Converts the ATAC cell-by-bin matrix into a gene activity matrix via `snapatac2.pp.make_gene_matrix` (hg38 GENCODE annotation).
10. **Modality tags** — GEX → `rna_unannotated`, ATAC activity → `atac_unannotated` (multiome only), reference → `rna_annotated`.
11. **Write outputs** — `~{input_id}_preprocessed_gex.h5ad` and `~{input_id}_preprocessed_ref.h5ad` always; `~{input_id}_preprocessed_atac_activity.h5ad` only in multiome mode.

#### 2. MultiomeLabelTransfer (GPU)

Loads the preprocessed h5ad files and performs **only** model training, label transfer, and output finalization. It imports individual functions (`run_multi_model`, `run_gex_only_model`, `transfer_labels`, `finalize_output`) from the container's `multiome_label_transfer.py` module — the script's `main()` function is **never called**, so no preprocessing is repeated.

1. **Load preprocessed data** — Reads the GEX and reference h5ad files produced by PreprocessFilter (plus ATAC activity in multiome mode). No filtering, reindexing, or conversion is performed.
2. **Train SCVI** — `run_multi_model()` (multiome) or `run_gex_only_model()` (GEX-only) concatenates the AnnData objects, filters to genes in ≥ 5 cells, selects 5,000 highly variable genes (Seurat v3, batch-aware), then trains an **SCVI** model (unsupervised VAE: 2 layers, 30 latent dimensions, negative-binomial likelihood, gene-batch dispersion, up to 500 epochs with early stopping).
3. **Train SCANVI** — The same function initializes **SCANVI** from the trained SCVI model and performs semi-supervised training using the reference cell type labels (`final_annotation`), propagating annotations to unlabeled query cells (up to 500 epochs, 100 samples per label).
4. **Predict labels** — `transfer_labels()` uses the trained SCANVI model to predict cell types (`C_scANVI`) for every cell, extracts the latent representation (`X_scANVI`), and computes a neighborhood graph and UMAP embedding.
5. **Propagate labels** — Copies the predicted `C_scANVI` labels from the concatenated object back into the original GEX (and, in multiome mode, ATAC) AnnData objects using the barcode-suffix index created by `ad.concat`.
6. **Write annotated matrices** — Saves `~{input_id}_gex_annotated_matrix.h5ad` (and `~{input_id}_atac_annotated_matrix.h5ad` in multiome mode).
7. **Finalize predictions** — `finalize_output()` adds placeholder metadata (biosample, donor, species, disease, organ, library prep, sex), renames `final_annotation` → `celltype`, and copies counts into the `.raw` layer for SCP ingest. Writes `~{input_id}_SCANVI_predictions.h5ad`.

### Workflow diagram

```text
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

In **GEX-only mode** (no ATAC input), the `_atac_activity` and `atac_annotated_matrix` branches are not produced; training uses `run_gex_only_model` on GEX + reference and labels are transferred to GEX only.

### Design rationale

The container script `multiome_label_transfer.py` bundles preprocessing **and** model training inside a single `main()` function. If the GPU task called `main()`, all preprocessing would run a second time on the GPU node — wasting expensive GPU hours on CPU-bound work that has already been completed.

Instead, the pipeline:

- Runs all preprocessing in **PreprocessFilter** on a CPU-only VM (no GPU cost).
- In **MultiomeLabelTransfer**, imports only the three model-training / label-transfer / finalization functions from the script, bypassing `main()` entirely.

This ensures zero duplication: every preprocessing step executes exactly once (on CPU), and the GPU node spends 100 % of its time on model training and inference.

### Outputs

All output filenames are prefixed with `~{input_id}_`.

| Output Variable Name | Filename | Output Type | Output Format |
| --- | --- | --- | --- |
| `scanvi_predictions_h5ad` | `<input_id>_SCANVI_predictions.h5ad` | Combined AnnData with SCANVI cell type predictions, UMAP, and metadata. | H5AD |
| `gex_annotated_h5ad` | `<input_id>_gex_annotated_matrix.h5ad` | Preprocessed GEX AnnData annotated with transferred cell type labels. | H5AD |
| `atac_annotated_h5ad` | `<input_id>_atac_annotated_matrix.h5ad` | ATAC gene-activity AnnData annotated with transferred cell type labels. Produced only in multiome mode (absent in GEX-only mode). | H5AD |
| `pipeline_version_out` | N/A | Version of the processing pipeline run on this data. | String |

## Runtime configuration

Both tasks use the same Docker image (pinned by digest). GPU and CUDA setup is handled entirely by the execution engine — the container does not configure the GPU environment itself.

#### Task 1 — `PreprocessFilter` (CPU-only)

| Attribute | Value |
| --- | --- |
| `docker` | `us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:635d4391d50cba9bd58f1fc41b10d8e1c61285a73bde75371815ce9a0db3430c` |
| `bootDiskSizeGb` | 20 |
| `disks` | `local-disk 1000 SSD` |
| `memory` | `120 GiB` |
| `cpu` | 32 |
| `maxRetries` | 1 |

#### Task 2 — `MultiomeLabelTransfer` (GPU)

| Attribute | Value |
| --- | --- |
| `docker` | `us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:635d4391d50cba9bd58f1fc41b10d8e1c61285a73bde75371815ce9a0db3430c` |
| `bootDiskSizeGb` | 20 |
| `disks` | `local-disk 500 SSD` |
| `memory` | `120 GiB` |
| `cpu` | 32 |
| `hardware_gpu_type` | `nvidia-tesla-t4` |
| `gpuCount` | 2 |
| `nvidia_driver_version` | `535.104.05` |
| `maxRetries` | 1 |

:::note GPU driver compatibility
Driver version `535.104.05` is compatible with CUDA 12.x and NVIDIA T4 GPUs and has been verified working on GCP / Terra with the `scvi-scanvi` container.
:::

## Docker image

The `scvi-scanvi` image is maintained in [warp-tools](https://github.com/broadinstitute/warp-tools/tree/develop/3rd-party-tools/scvi-scanvi). Key libraries: scvi-tools 1.2, snapatac2 2.7, scanpy, anndata.

## Versioning

All scANVI pipeline releases are documented in the [scANVI changelog](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.changelog.md).

## Citing the scANVI Pipeline

When citing WARP, please use the following:

Kylee Degatano, Aseel Awdeh, Robert Sidney Cox III, Wes Dingman, George Grant, Farzaneh Khajouei, Elizabeth Kiernan, Kishori Konwar, Kaylee L Mathews, Kevin Palis, Nikelle Petrillo, Geraldine Van der Auwera, Chengchen (Rex) Wang, Jessica Way. "Warp Analysis Research Pipelines: Cloud-optimized workflows for biological data processing and reproducible analysis." _Bioinformatics_, 2025; [https://doi.org/10.1093/bioinformatics/btaf494](https://doi.org/10.1093/bioinformatics/btaf494)

Please also cite the underlying scvi-tools models:

- Lopez, R., Regier, J., Cole, M.B., Jordan, M.I., Yosef, N. "Deep generative modeling for single-cell transcriptomics." _Nature Methods_ 15, 1053–1058 (2018).
- Xu, C., Lopez, R., Mehlman, E., Regier, J., Jordan, M.I., Yosef, N. "Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models." _Molecular Systems Biology_ 17, e9620 (2021).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.
