---
sidebar_position: 1
slug: /Pipelines/scANVI_Pipeline/README
---

# scANVI Overview

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| See [changelog](https://github.com/broadinstitute/warp/blob/develop/pipelines/wdl/scanvi/scANVI.changelog.md) for version information. | See changelog | WARP Pipelines | Please [file an issue in WARP](https://github.com/broadinstitute/warp/issues)  |

## Introduction to the scANVI workflow

The scANVI pipeline is a cloud-optimized WDL workflow that performs **cell type label transfer** using [scVI](https://docs.scvi-tools.org/en/1.4.1/user_guide/models/scvi.html) and [scANVI](https://docs.scvi-tools.org/en/1.4.1/user_guide/models/scanvi.html) (single-cell ANnotation using Variational Inference) deep generative models. It integrates single-cell RNA-seq (GEX) and, optionally, ATAC-seq data with an annotated reference dataset to transfer cell type labels via semi-supervised learning.

ATAC is optional. When no ATAC h5ads are supplied, the pipeline runs in **GEX-only mode** and trains/annotates from the reference atlas using gene expression and reference data alone, without using ATAC. When ATAC is supplied it runs in multiome mode (GEX + ATAC + reference).

One run labels **one or more query files** against a single annotated reference. The reference is prepared once; each query is then preprocessed and labelled independently, in parallel, and gets its own outputs.

The pipeline is split into three tasks: **PreprocessReference** (CPU-only, once per run), **PreprocessFilter** (CPU-only, once per query) and **MultiomeLabelTransfer** (GPU, once per query), so that expensive GPU time is reserved exclusively for model training and inference. All data loading, quality filtering, and (in multiome mode) barcode alignment and gene-activity-matrix conversion happen on a CPU node first; the GPU node receives ready-to-train h5ad files and never re-runs preprocessing.

### How the label transfer works: SCVI + SCANVI

scANVI annotates in two stages — an **unsupervised** representation-learning step followed by a **semi-supervised** annotation step — so it can exploit *all* the cells (the large unlabeled query alongside the annotated reference), not just the labeled ones.

1. **SCVI — unsupervised embedding.** [SCVI](https://docs.scvi-tools.org/en/1.4.1/user_guide/models/scvi.html) is a variational autoencoder trained on the **concatenation** of the query and reference cells (GEX, plus the ATAC gene-activity matrix in multiome mode) **without using any cell-type labels**. It learns a low-dimensional latent representation that captures biological variation while correcting batch effects — across donors, and between the query and the reference (and across modalities). Being unsupervised, it uses every cell, so the embedding is shaped by the full dataset rather than only the annotated subset.

2. **SCANVI — semi-supervised annotation.** [SCANVI](https://docs.scvi-tools.org/en/1.4.1/user_guide/models/scanvi.html) is initialized **from the trained SCVI model** and adds a cell-type classifier on top of that latent space. It is trained **semi-supervised**: the reference cells provide their known labels (the supervision signal) while the unlabeled query cells (tagged `Unknown`) continue to shape the latent space (the unsupervised signal). The result is a label-aware embedding in which query cells sit near reference cells of the same type, and SCANVI then predicts a cell type for every query cell (`C_scANVI`).

This unsupervised → semi-supervised design is more robust than training a supervised classifier on the reference alone: SCVI first integrates query and reference into a common, batch-corrected space using all available cells, and SCANVI only has to learn the annotation boundaries **within** that shared space. Labels are transferred by propagating each query cell's `C_scANVI` prediction back onto the original matrices.

#### ATAC gene-activity conversion (Multiome mode)

Because SCVI requires all integrated datasets to share the exact same feature space, the raw ATAC cell-by-bin matrix (produced by upstream pipelines) cannot be concatenated directly with the cell-by-gene RNA matrices. The pipeline solves this by projecting the ATAC signal into the transcriptomic feature space:

1. **Bypassing peak calling.** Rather than relying on distal peaks (which are difficult to definitively map to target genes), the pipeline counts the raw Tn5 insertion fragments that overlap each known gene's coordinates.
2. **Gene activity scoring.** Using a genome annotation (e.g., hg38 GENCODE), `snapatac2.pp.make_gene_matrix` defines a regulatory domain for each gene—typically the gene body plus a promoter window immediately upstream of the transcription start site. The accessibility metric for that gene is the count of fragments falling within that window; fragments falling in distant intergenic regions are ignored.
3. **Concatenation.** This produces a cell-by-gene "activity" matrix that mathematically resembles an RNA count matrix. The GEX query, the ATAC activity matrix, and the reference can then be seamlessly concatenated, reduced to the top highly variable genes, and mapped into a single SCVI latent space where technical modality differences are handled mathematically as a batch effect.

#### How the models are trained

Both models are trained by **minibatch stochastic gradient descent**. The full concatenated dataset (query + reference) is prepared and held in CPU/host memory, but each training step streams only a small **minibatch** of cells — scvi-tools' default `batch_size` is **128 cells** — to the GPU, runs the forward/backward pass, updates the model weights, and releases that minibatch before fetching the next. Because the GPU only ever holds one minibatch at a time, its memory footprint is set by *minibatch size × number of genes* and is essentially **independent of the total number of cells** in the dataset. This is what lets a large query be annotated on a single modest GPU while the full dataset lives in host RAM: the CPU node assembles and holds the data, and the GPU processes it 128 cells at a time.

**Two meanings of "batch."** The `batch_size` above is the SGD minibatch — an optimization detail of how the data is fed to the GPU. It is distinct from the **batch covariate**, the experimental grouping (e.g., donor or sequencing library) that SCVI/SCANVI explicitly model in order to correct for it as technical variation. The models integrate *across* batch-covariate groups to remove batch effects, and they do so by reading the data in `batch_size`-cell minibatches — the same word, two unrelated concepts.

:::tip Want to use scANVI for your publication?
The pipeline is designed to consume the outputs of the [Multiome](../Multiome_Pipeline/README.md) and [PeakCalling](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/peak_calling) WARP pipelines. Identify it in your methods section with its SciCrunch resource identifier, *scANVI Pipeline (RRID:SCR_028705)*, and cite the pipeline using the WARP citation in the [Citing](#citing-the-scanvi-pipeline) section below.
:::

## Quickstart table

The following table provides a quick glance at the scANVI pipeline features:

| Pipeline features | Description | Source |
| --- | --- | --- |
| Assay type | 10x single-cell / single-nucleus Multiome (GEX + ATAC), or GEX-only | [10x Genomics](https://www.10xgenomics.com) |
| Overall workflow | CPU reference preparation (once per run), then per-query CPU preprocessing and GPU SCVI/SCANVI label transfer | Code available on [GitHub](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.wdl) |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Models | SCVI (unsupervised VAE) + SCANVI (semi-supervised classifier) | [scvi-tools 1.5.1](https://docs.scvi-tools.org/) |
| ATAC gene-activity conversion | Cell-by-bin matrix → gene activity matrix (hg38 GENCODE) | [snapatac2 2.10](https://kzhang.org/SnapATAC2/) |
| Data input format | AnnData h5ad files: one annotated reference and one or more GEX queries (required), plus one ATAC cell-by-bin file per query (optional) | [AnnData](https://anndata.readthedocs.io/) |
| Data output format | Annotated h5ad files with predicted cell types and UMAP | [AnnData](https://anndata.readthedocs.io/) |

## Set-up

### scANVI installation

To download the latest scANVI release, see the release tags prefixed with "scANVI" on the WARP [releases page](https://github.com/broadinstitute/warp/releases). All scANVI pipeline releases are documented in the [scANVI changelog](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.changelog.md).

scANVI can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant, flexible workflow management system that supports multiple computing platforms. The workflow can also be run in [Terra](https://app.terra.bio).

### Inputs

scANVI takes its queries as arrays paired by position: `input_ids[i]`, `gex_h5ads[i]` and, in multiome mode, `atac_h5ads[i]` describe query *i*, and every query is labelled against the same `ref_h5ad`. A single query is a one-element array.

Example input JSON files are available in the [`example_inputs`](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/scanvi/example_inputs) folder.

#### Workflow inputs

| Input | Type | Description | Default |
| --- | --- | --- | --- |
| `input_ids` | Array[String] | One identifier per query, in the same order as `gex_h5ads`; `input_ids[i]` is prepended to every output filename of query *i*, so ids must be unique. | — (required) |
| `gex_h5ads` | Array[File] | Gene expression AnnData h5ad files from Multiome / Optimus, one per query. | — (required) |
| `atac_h5ads` | Array[File] | **Optional** ATAC cell-by-bin AnnData h5ad files from Multiome / PeakCalling, one per `gex_h5ads` entry and in the same order; empty runs every query in GEX-only mode. | `[]` |
| `ref_h5ad` | File | Annotated reference AnnData h5ad file, shared by every query. | — (required) |
| `max_epochs` | Int? | Optional cap on SCVI/SCANVI training epochs, applied in both multiome and GEX-only modes. When unset, the container default (500) is used. | — |
| `batch_size` | Int | SCVI/SCANVI minibatch size (the SGD minibatch described above). Lower it to fit a high-cardinality reference on a small GPU (activation memory scales with `batch_size` × number of labels); raise it on a large-VRAM cloud GPU. | `128` |
| `seed` | Int? | Optional random seed for SCVI/SCANVI training. Set it to make a run reproducible (the test wrappers pin one); unset keeps training nondeterministic. | — |
| `ref_label_column` | String? | Reference `obs` column to use as the cell-type label. When unset, defaults to `subclass` for AIT references and `final_annotation` otherwise. | — |
| `ref_batch_column` | String? | Reference `obs` column to use as the batch. When unset, defaults to `donor_id` for AIT references and `batch` otherwise. | — |
| `genome` | String | Genome for the ATAC cell-by-bin → gene-activity conversion (multiome only): `hg38` (default), `mm10`, or `mm39`. | `"hg38"` |
| `scanvi_model` | File? | **Optional** pre-trained SCANVI model (`.tar.gz` of a saved model directory, no bundled AnnData — the format emitted in `scanvi_models_out`). When provided, the label-transfer task **loads it and predicts every query with it, skipping SCVI/SCANVI training** (auto-detected, no flag). Valid when the model matches the incoming data/reference. Without it, each query trains its own model. | — |
| `output_max_probability` | Boolean | When `true`, add a `max_probability` obs column (per-cell maximum SCANVI posterior — the assigned label's confidence) to every output h5ad. | `false` |
| `gpu_count` | Int | GPUs requested for each query's label-transfer task. Set `0` for a CPU-only prediction run (e.g. a supplied-model run on small data). | `2` |
| `mem_size` | Int | Memory (GiB) for each query's label-transfer task. | `120` |
| `nthreads` | Int | CPUs for each query's label-transfer task. | `32` |
| `disk_size` | Int | Disk (GB) for each query's label-transfer task. | `500` |

:::note One query or many

- Pairing is by position, and the run fails before any query is processed if `gex_h5ads` is empty, `input_ids` differs in length, or `atac_h5ads` is neither empty nor the same length. ATAC is all or none per run.
- In an inputs JSON a single query still needs brackets: `"scANVI.input_ids": ["10k_PBMC"]`, `"scANVI.gex_h5ads": ["gs://…"]`. A Terra data-table column such as `this.gex_h5ad` maps to an array input directly, because Terra wraps the single value.
- Each query's outputs are named after its input id, so ids must be unique.
- Without `scanvi_model` each query trains its own model, so cost grows with the number of queries. With it, every query is predicted by the same frozen model.

:::

:::note GEX-only mode
The ATAC input is optional. Leave `atac_h5ads` empty (the default) and every query runs in GEX-only mode: SCVI/SCANVI are trained on GEX + reference via `run_gex_only_model`, labels are transferred to GEX only, and `atac_annotated_h5ads` is empty. See [`example_inputs/scANVI.gex_only.json`](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/scanvi/example_inputs).
:::

:::note Reusing a trained model (skip training)
Every query emits its SCANVI model in `scanvi_models_out` — a small `.tar.gz` (no bundled data). Pass one back as the `scanvi_model` input on a later run and every query is predicted with it: the label-transfer task **loads it and predicts instead of training SCVI/SCANVI**, which is valid when the model matches the incoming data/reference. For a cheap CPU-only prediction run, also set `gpu_count = 0` (and smaller `mem_size`/`nthreads`/`disk_size`). Training once and then predicting many queries on CPU this way is far cheaper than training per query (see [Example run costs](#example-run-costs)). Applying a reference model to genuinely novel query data (scArches-style query mapping) is planned future work.
:::

#### Running on Terra: per sample or per set

Terra can feed the query arrays from a data table in two ways. The column names below (`sample_id`, `gex_h5ad`, `atac_h5ad`, `ref_h5ad`) are examples; use your own table's names.

- **Per row** (root entity `sample`): set `input_ids` = `this.sample_id`, `gex_h5ads` = `this.gex_h5ad`, `atac_h5ads` = `this.atac_h5ad` and `ref_h5ad` = `workspace.ref_h5ad`. Each selected row runs as its own workflow, and Terra wraps each single value into a one-element array. A row with no ATAC value runs GEX-only. Outputs land on the row as one-element lists. Rows launched together each prepare the reference themselves.
- **Per set** (root entity `sample_set`): set `input_ids` = `this.samples.sample_id`, `gex_h5ads` = `this.samples.gex_h5ad` and `atac_h5ads` = `this.samples.atac_h5ad` (leave it unset for a GEX-only set), with `ref_h5ad` as above. One workflow runs per set and the reference is prepared once. A supplied model is used for every sample. Outputs land on the set row in sample order. A set that mixes multiome and GEX-only samples fails the pairing check, so split it.

Without a data table, list the queries in the inputs JSON. For two multiome queries (placeholder paths):

```json
{
  "scANVI.input_ids": ["sample_a", "sample_b"],
  "scANVI.gex_h5ads": [
    "gs://your-bucket/sample_a_gex.h5ad",
    "gs://your-bucket/sample_b_gex.h5ad"
  ],
  "scANVI.atac_h5ads": [
    "gs://your-bucket/sample_a_atac.h5ad",
    "gs://your-bucket/sample_b_atac.h5ad"
  ],
  "scANVI.ref_h5ad": "gs://your-bucket/reference.h5ad"
}
```

#### Reference requirements

The reference h5ad must provide cell-type labels and a batch. Two reference forms are supported:

- **PBMC-style** (default): counts in `.X`, labels in `obs['final_annotation']`, batch in `obs['batch']`.
- **AIT-schema** (Allen Institute Taxonomy): auto-detected via `uns['schema_version']` + `uns['hierarchy']`. These typically have **no `.X`** (counts are materialized from `.raw`, keeping the gene symbols from `.var`), a cell-type **hierarchy** (e.g. `class`/`subclass`/`cluster_id`), and a batch column such as `donor_id`. By default the label is taken from `subclass` and the batch from `donor_id`; override with `ref_label_column` / `ref_batch_column`. For a mouse AIT reference used in multiome mode, set `genome` to `mm10`/`mm39` so the ATAC gene-activity conversion uses the correct genome.

The query datasets (GEX, and ATAC when present) do not need pre-existing annotations — placeholder `Unknown` labels are added automatically before training.

## scANVI tasks and tools

The [scANVI workflow](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.wdl) defines three tasks inline, plus `MultiomeLabelTransferCpu`, the CPU-only twin of MultiomeLabelTransfer. All use the same Docker image, and only MultiomeLabelTransfer is allocated GPUs. PreprocessReference runs once per run; the other two run once per query.

| Task name | Tool | Software | Description |
| --- | --- | --- | --- |
| [PreprocessReference](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.wdl) | scanpy | Python | Prepares the reference once per run: materializes counts from `.raw` for AIT references, maps the chosen label/batch columns to `final_annotation`/`batch`, and tags the reference cells with modality `rna_annotated`. |
| [PreprocessFilter](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.wdl) | scanpy, snapatac2 | Python | Loads one query's GEX h5ad (plus its ATAC h5ad in multiome mode), patches missing columns, and filters GEX to STARsolo cell calls. In multiome mode it also intersects barcodes between GEX and ATAC and converts the ATAC cell-by-bin matrix into a gene activity matrix; in GEX-only mode these ATAC steps are skipped. |
| [MultiomeLabelTransfer](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.wdl) | scvi-tools | Python (GPU) | Trains an SCVI model on the preprocessed AnnData objects (GEX + ATAC + reference, or GEX + reference in GEX-only mode), then trains an SCANVI classifier from the SCVI model to transfer reference cell-type labels onto the unlabeled query cells. Computes UMAP and writes annotated outputs. |

Overall, the scANVI workflow:

1. [Prepares the annotated reference once per run (CPU).](#1-preprocessreference-cpu-only)
2. [Preprocesses and filters each query's h5ad files (CPU).](#2-preprocessfilter-cpu-only)
3. [Trains SCVI / SCANVI models and transfers labels, once per query (GPU).](#3-multiomelabeltransfer-gpu)

#### 1. PreprocessReference (CPU-only)

Prepares the annotated reference on a CPU-only node. It runs once per run, outside the per-query scatter, and every query shares its output. Steps:

1. **Load the reference** — Reads `ref_h5ad` (scanpy).
2. **Adapt the reference** — Detects an AIT-schema reference from `uns['schema_version']` plus `uns['hierarchy']`; materializes counts from `.raw` (keeping the `.var` gene symbols) when `.X` is absent; and maps the label and batch columns, from `ref_label_column` / `ref_batch_column` or the defaults (see [Reference requirements](#reference-requirements)), to `final_annotation` and `batch`.
3. **Modality tag** — Tags the reference cells `rna_annotated`.
4. **Write output** — `preprocessed_ref.h5ad`, which every query's label-transfer task reads.

#### 2. PreprocessFilter (CPU-only)

Loads and preprocesses one query's h5ad files on a CPU-only node. It runs once per query, inside the per-query scatter. No GPU is allocated for this task. Steps 5, 6, and 9 (the ATAC-specific steps) run **only in multiome mode**; in GEX-only mode they are skipped and no ATAC activity file is written. Steps:

1. **Load datasets** — Reads the query's GEX (scanpy), plus its ATAC cell-by-bin (snapatac2) in multiome mode.
2. **Patch missing columns** — Adds `star_IsCell = True` to GEX and `gex_barcodes` (from index) to ATAC if absent, ensuring compatibility across upstream pipelines.
3. **Filter GEX** — Retains STARsolo cell calls (`star_IsCell == True`), then removes genes and cells with fewer than 3 counts.
4. **Prepare GEX** — Sets `batch` label; copies counts into a `counts` layer.
5. **Reindex ATAC** *(multiome only)* — Sets ATAC obs index to `gex_barcodes` so barcodes align with GEX.
6. **Shared barcode filtering** *(multiome only)* — Intersects GEX and ATAC barcodes; subsets both to matched cells. In GEX-only mode all filtered GEX cells are retained.
7. **Batch labels** — GEX → `pd-multiome_sci_gex` (and, in multiome mode, ATAC → `pd-multiome_sci_atac`).
8. **Placeholder annotations** — Adds `final_annotation = "Unknown"` to query datasets.
9. **Gene activity matrix** *(multiome only)* — Converts the ATAC cell-by-bin matrix into a gene activity matrix via `snapatac2.pp.make_gene_matrix` (hg38 GENCODE annotation).
10. **Modality tags** — GEX → `rna_unannotated`, ATAC activity → `atac_unannotated` (multiome only).
11. **Write outputs** — `~{input_id}_preprocessed_gex.h5ad` always; `~{input_id}_preprocessed_atac_activity.h5ad` only in multiome mode.

#### 3. MultiomeLabelTransfer (GPU)

Runs once per query, inside the per-query scatter. Loads the preprocessed h5ad files and performs **only** model training, label transfer, and output finalization. It imports individual functions (`run_multi_model`, `run_gex_only_model`, `transfer_labels`, `finalize_output`) from the container's `multiome_label_transfer.py` module — the script's `main()` function is **never called**, so no preprocessing is repeated.

1. **Load preprocessed data** — Reads the query's preprocessed GEX (plus ATAC activity in multiome mode) from PreprocessFilter and the preprocessed reference from PreprocessReference. No filtering, reindexing, or conversion is performed.
2. **Train SCVI** — `run_multi_model()` (multiome) or `run_gex_only_model()` (GEX-only) concatenates the AnnData objects, filters to genes in ≥ 5 cells, selects 5,000 highly variable genes (Seurat v3, batch-aware), then trains an **SCVI** model (unsupervised VAE: 2 layers, 30 latent dimensions, negative-binomial likelihood, gene-batch dispersion, up to 500 epochs with early stopping).
3. **Train SCANVI** — The same function initializes **SCANVI** from the trained SCVI model and performs semi-supervised training using the reference cell type labels (`final_annotation`), propagating annotations to unlabeled query cells (up to 500 epochs, 100 samples per label).
4. **Predict labels** — `transfer_labels()` uses the trained SCANVI model to predict cell types (`C_scANVI`) for every cell, extracts the latent representation (`X_scANVI`), and computes a neighborhood graph and UMAP embedding.
5. **Propagate labels** — Copies the predicted `C_scANVI` labels from the concatenated object back into the original GEX (and, in multiome mode, ATAC) AnnData objects using the barcode-suffix index created by `ad.concat`.
6. **Write annotated matrices** — Saves `~{input_id}_gex_annotated_matrix.h5ad` (and `~{input_id}_atac_annotated_matrix.h5ad` in multiome mode).
7. **Finalize predictions** — `finalize_output()` renames `final_annotation` → `celltype` and copies counts into the `.raw` layer for downstream ingest. No sample metadata is added. Writes `~{input_id}_SCANVI_predictions.h5ad`.

### Workflow diagram

```text
        ref_h5ad
            │
            ▼
┌───────────────────────┐
│  PreprocessReference  │  CPU-only, once per run
│ (adapt the reference) │
└───────────┬───────────┘
            │ preprocessed_ref.h5ad (shared by every query)
            │
┌───────────┼─ for each query i (scatter, in parallel) ─────────┐
│           │                                                   │
│           │  gex_h5ads[i] (+ atac_h5ads[i])                   │
│           │                 │                                 │
│           │                 ▼                                 │
│           │    ┌─────────────────────────┐                    │
│           │    │    PreprocessFilter     │  CPU-only          │
│           │    │  (load, filter, align,  │                    │
│           │    │  gene activity matrix)  │                    │
│           │    └────────────┬────────────┘                    │
│           │                 │ preprocessed GEX                │
│           │                 │ (+ ATAC activity)               │
│           │                 │                                 │
│           ▼                 ▼                                 │
│     ┌───────────────────────────────────────────────┐         │
│     │             MultiomeLabelTransfer             │  GPU    │
│     │ (MultiomeLabelTransferCpu when gpu_count = 0) │         │
│     └───────────────────────┬───────────────────────┘         │
│                             │                                 │
└─────────────────────────────┼─────────────────────────────────┘
                              │
                              ▼
           input_ids[i]_SCANVI_predictions.h5ad
           input_ids[i]_gex_annotated_matrix.h5ad
           input_ids[i]_atac_annotated_matrix.h5ad   (multiome only)
           input_ids[i]_scanvi_model.tar.gz
```

In **GEX-only mode** (empty `atac_h5ads`), the `_atac_activity` and `atac_annotated_matrix` branches are not produced; training uses `run_gex_only_model` on GEX + reference and labels are transferred to GEX only.

### Design rationale

The container script `multiome_label_transfer.py` bundles preprocessing **and** model training inside a single `main()` function. If the GPU task called `main()`, all preprocessing would run a second time on the GPU node — wasting expensive GPU hours on CPU-bound work that has already been completed.

Instead, the pipeline:

- Runs all preprocessing in **PreprocessReference** and **PreprocessFilter** on CPU-only VMs (no GPU cost).
- In **MultiomeLabelTransfer**, imports only the three model-training / label-transfer / finalization functions from the script, bypassing `main()` entirely.
- Prepares the shared reference once, in **PreprocessReference** outside the per-query scatter, so a multi-query run loads, adapts and writes it once instead of once per query, and call caching can reuse it.

This ensures zero duplication: every preprocessing step executes exactly once (on CPU), and the GPU node spends 100 % of its time on model training and inference.

### Outputs

Every output except `pipeline_version_out` is an array with one entry per query, in the order of `input_ids`; each file name starts with that query's input id.

| Output Variable Name | Filename | Output Type | Output Format |
| --- | --- | --- | --- |
| `scanvi_predictions_h5ads` | `<input_id>_SCANVI_predictions.h5ad` | Combined AnnData (query + reference cells) with SCANVI cell type predictions (`celltype`) and UMAP; counts in `.raw`. | H5AD |
| `gex_annotated_h5ads` | `<input_id>_gex_annotated_matrix.h5ad` | Preprocessed GEX AnnData annotated with transferred cell type labels. | H5AD |
| `atac_annotated_h5ads` | `<input_id>_atac_annotated_matrix.h5ad` | ATAC gene-activity AnnData annotated with transferred cell type labels. Produced only in multiome mode; empty in GEX-only mode. | H5AD |
| `scanvi_models_out` | `<input_id>_scanvi_model.tar.gz` | The query's trained (or loaded) SCANVI model, saved without bundled AnnData (~tens of MB). Can be passed back as `scanvi_model` on a later run to skip training. | TAR.GZ |
| `pipeline_version_out` | N/A | Version of the processing pipeline run on this data. | String |

:::note No sample metadata
The pipeline adds no sample-level metadata (donor, species, organ, disease, library preparation, sex) to any output. The predictions file retains only input `obs` columns shared by all concatenated datasets, plus pipeline-generated annotation columns; columns present in only one input, such as a reference-only `donor_id`, are not retained there. The annotated matrices retain their preprocessed modality-specific `obs` columns and transferred labels. Sample metadata is joined downstream from the data archive's sample manifests, not by this pipeline.
:::

## Runtime configuration

All tasks use the same Docker image (pinned by digest). GPU and CUDA setup is handled entirely by the execution engine — the container does not configure the GPU environment itself.

#### Task 1 — `PreprocessReference` (CPU-only, once per run)

| Attribute | Value |
| --- | --- |
| `docker` | `us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:077d776c1957029b32ef79f12b6dac02e43d5673b70e915d1d690d95612a9831` |
| `bootDiskSizeGb` | 20 |
| `disks` | `local-disk 1000 SSD` |
| `memory` | `120 GiB` |
| `cpu` | 32 |
| `maxRetries` | 1 |

#### Task 2 — `PreprocessFilter` (CPU-only, once per query)

| Attribute | Value |
| --- | --- |
| `docker` | `us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:077d776c1957029b32ef79f12b6dac02e43d5673b70e915d1d690d95612a9831` |
| `bootDiskSizeGb` | 20 |
| `disks` | `local-disk 1000 SSD` |
| `memory` | `120 GiB` |
| `cpu` | 32 |
| `maxRetries` | 1 |

#### Task 3 — `MultiomeLabelTransfer` (GPU, once per query)

| Attribute | Value |
| --- | --- |
| `docker` | `us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:077d776c1957029b32ef79f12b6dac02e43d5673b70e915d1d690d95612a9831` |
| `bootDiskSizeGb` | 20 |
| `disks` | `local-disk 500 SSD` |
| `memory` | `120 GiB` |
| `cpu` | 32 |
| `gpuType` | `nvidia-tesla-t4` |
| `gpuCount` | `gpu_count` (default 2) |
| `nvidiaDriverVersion` | `535.104.05` |
| `maxRetries` | 1 |

:::note GPU driver compatibility
Driver version `535.104.05` is compatible with CUDA 12.x and NVIDIA T4 GPUs and has been verified working on GCP / Terra with the `scvi-scanvi` container.
:::

:::note CPU-only variant
When `gpu_count = 0` the workflow routes to `MultiomeLabelTransferCpu` — an identical task with **no** GPU runtime attributes (Cromwell rejects `gpuCount = 0` and can't conditionally omit the GPU keys from one task). Both tasks run the same container command (`label_transfer_from_preprocessed.py`); scvi-tools auto-detects the accelerator, so it uses the GPU when present and CPU otherwise.
:::

## Example run costs

The table below reports representative Terra costs for scANVI runs across three references — mouse hippocampus, human neocortex, and PBMC — in both training and CPU-inference modes. Costs are from single runs on GCP `us-central1` (training on **2× NVIDIA T4** GPUs; inference CPU-only via `gpu_count = 0`) and will vary with region, machine type, preemption, and dataset. "Cells labelled" is the number of query cells annotated (after PreprocessFilter; for multiome, after the GEX↔ATAC shared-barcode intersection); "Reference" is the annotated `ref_h5ad`.

Each row is a single-query run. A multi-query run prepares the reference once, then runs PreprocessFilter and label transfer once per query. Without `scanvi_model` each query also trains its own model, so a multi-query training run costs roughly the train cost for that reference times the number of queries. Predicting every query with one saved model on CPU (`scanvi_model` plus `gpu_count = 0`) is far cheaper.

| Dataset (reference) | Mode | Compute | Cost | Cells labelled | Cost / cell | Cost / 1k cells | Query input | Cost / GB (query) | Reference (cells / GB) |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Mouse hippocampus (AIT) | train | 2× T4 GPU | $26.34 | 24,078 | $0.001094 | $1.09 | 0.53 GB | $49.85 | 250,734 / 12.87 |
| Mouse hippocampus (AIT) | inference | CPU | $0.56 | 24,078 | $0.0000233 | $0.02 | 1.69 GB | $0.33 | 250,734 / 12.87 |
| Human neocortex (AIT) | train | 2× T4 GPU | $2.51 | 6,018 | $0.000417 | $0.42 | 0.32 GB | $7.98 | 47,432 / 4.97 |
| Human neocortex (AIT) | inference | CPU | $0.33 | 6,018 | $0.0000548 | $0.05 | 0.32 GB | $1.05 | 47,432 / 4.97 |
| 10k PBMC (multiome) | train | 2× T4 GPU | $3.64 | 911 | $0.003996 | $4.00 | 0.04 GB GEX + 0.96 GB ATAC | $100.14 † | 33,506 / 2.06 |
| 10k PBMC (GEX-only) | train | 2× T4 GPU | $3.38 | 934 | $0.003619 | $3.62 | 0.04 GB | $92.98 | 33,506 / 2.06 |

† GEX file only; $3.66/GB when counting GEX + ATAC input together.

Takeaways:

- **Inference is far cheaper than training.** Loading a saved model and predicting on CPU (`scanvi_model` supplied, `gpu_count = 0`) cost 8–47× less than the equivalent training run.
- **Training cost tracks the reference, not the query.** The largest query (1.69 GB, mouse inference) was among the cheapest runs; training cost scales with reference size (mouse 250,734-cell / 12.87 GB → $26.34 vs human 47,432 / 4.97 GB → $2.51).
- **Per-cell cost rises as the query shrinks**, because the fixed reference-load / training cost is amortized over fewer labelled cells (e.g. PBMC's ~900-cell queries are the highest cost/cell).
- **Multiome and GEX-only cost about the same** ($3.64 vs $3.38) — the ATAC branch adds little.

## Docker image

The `scvi-scanvi` image is maintained in [warp-tools](https://github.com/broadinstitute/warp-tools/tree/develop/3rd-party-tools/scvi-scanvi). Key libraries: scvi-tools 1.5.1, snapatac2 2.10, scanpy, anndata.

## Versioning

All scANVI pipeline releases are documented in the [scANVI changelog](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.changelog.md).

## Citing the scANVI Pipeline

If you use the scANVI Pipeline in your research, please identify the pipeline in your methods section using the [scANVI SciCrunch resource identifier](https://rrid.site/resolver/SCR_028705).

- Ex: *scANVI Pipeline (RRID:SCR_028705)*

When citing WARP, please use the following:

Kylee Degatano, Aseel Awdeh, Robert Sidney Cox III, Wes Dingman, George Grant, Farzaneh Khajouei, Elizabeth Kiernan, Kishori Konwar, Kaylee L Mathews, Kevin Palis, Nikelle Petrillo, Geraldine Van der Auwera, Chengchen (Rex) Wang, Jessica Way. "Warp Analysis Research Pipelines: Cloud-optimized workflows for biological data processing and reproducible analysis." _Bioinformatics_, 2025; [https://doi.org/10.1093/bioinformatics/btaf494](https://doi.org/10.1093/bioinformatics/btaf494)

Please also cite the underlying scvi-tools models:

- Lopez, R., Regier, J., Cole, M.B., Jordan, M.I., Yosef, N. "Deep generative modeling for single-cell transcriptomics." _Nature Methods_ 15, 1053–1058 (2018).
- Xu, C., Lopez, R., Mehlman, E., Regier, J., Jordan, M.I., Yosef, N. "Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models." _Molecular Systems Biology_ 17, e9620 (2021).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.
