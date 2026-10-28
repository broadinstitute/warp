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

One run labels **one or more samples**. Each sample is labelled independently and gets its own outputs, and every sample of a run goes through **one machine, one after another**, so a run of many samples starts one machine instead of one per sample. The model comes from one of three places:

- **No model supplied:** each sample trains its own SCVI/SCANVI model on the annotated reference plus that sample.
- **A model this pipeline trained** (`scanvi_models_out` of an earlier run): it is loaded and predicts every sample as is.
- **A model trained on the reference alone**, for example by [warp-scanvi-training](https://github.com/rsc3/warp-scanvi-training): each sample is mapped onto it with scArches (query mapping), a short training on the sample alone.

With a supplied model the reference is not needed. The annotated matrices are written ready for [Single Cell Portal](https://singlecell.broadinstitute.org) (SCP): counts in `.raw`, normalized expression in `.X`, a 2-D UMAP and the scANVI latent. Every output records which sample, input file and model produced it.

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

scANVI takes its samples as arrays paired by position: `input_ids[i]`, `gex_h5ads[i]` and, in multiome mode, `atac_h5ads[i]` describe sample *i*. A single sample is a one-element array. Give either `scanvi_model` (to label with a saved model) or `ref_h5ad` (to train); the run fails at once with neither.

Example input JSON files are available in the [`example_inputs`](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/scanvi/example_inputs) folder.

#### Workflow inputs

| Input | Type | Description | Default |
| --- | --- | --- | --- |
| `input_ids` | Array[String] | One identifier per sample, in the same order as `gex_h5ads`. `input_ids[i]` starts every output file name of sample *i* and is recorded in each of its outputs, so ids must be unique and cannot contain `/`. | — (required) |
| `gex_h5ads` | Array[File] | Gene expression AnnData h5ad files from Multiome / Optimus, one per sample. | — (required) |
| `atac_h5ads` | Array[File] | **Optional** ATAC cell-by-bin AnnData h5ad files from Multiome / PeakCalling, one per `gex_h5ads` entry and in the same order; empty runs every sample in GEX-only mode. | `[]` |
| `scanvi_model` | File? | **Optional** saved SCANVI model (`.tar.gz` of a saved model directory, no bundled AnnData), used for every sample instead of training. A model this pipeline trained loads as is; a model trained on the reference alone gets each sample mapped onto it (see [Using a saved model](#using-a-saved-model)). | — |
| `ref_h5ad` | File? | Annotated reference AnnData h5ad file. Needed only to train, that is, when `scanvi_model` is not given; ignored otherwise. | — |
| `max_epochs` | Int? | Optional cap on SCVI/SCANVI training epochs, applied in both multiome and GEX-only modes. When unset, the container default (500) is used. | — |
| `query_max_epochs` | Int? | Epochs when mapping a sample onto a reference-only model. When unset, the container default (100) is used. | — |
| `batch_size` | Int | SCVI/SCANVI minibatch size (the SGD minibatch described above). Lower it to fit a high-cardinality reference on a small GPU (activation memory scales with `batch_size` × number of labels); raise it on a large-VRAM cloud GPU. | `128` |
| `seed` | Int? | Optional random seed for SCVI/SCANVI training. Set it to make a run reproducible (the test wrappers pin one); unset keeps training nondeterministic. | — |
| `ref_label_column` | String? | When training, the reference `obs` column to use as the cell-type label. When unset, defaults to `subclass` for AIT references and `final_annotation` otherwise. | — |
| `ref_batch_column` | String? | When training, the reference `obs` column to use as the batch. When unset, defaults to `donor_id` for AIT references and `batch` otherwise. | — |
| `genome` | String | Genome for the ATAC cell-by-bin → gene-activity conversion (multiome only): `hg38` (default), `mm10`, or `mm39`. | `"hg38"` |
| `output_max_probability` | Boolean | When `true`, add a `max_probability` obs column (per-cell maximum SCANVI posterior — the assigned label's confidence) to every output h5ad. | `false` |
| `gpu_count` | Int | GPUs for the label-transfer machine. Set `0` for a CPU-only run, which suits a supplied model. | `2` |
| `mem_size` | Int | Memory (GiB) for the label-transfer machine. Samples run one at a time, so it only has to fit the largest sample (plus the reference when training). | `104` |
| `nthreads` | Int | CPUs for the label-transfer machine. | `16` |
| `disk_size` | Int? | Disk (GB) for the label-transfer machine. When unset, it is sized from the inputs: 4 × the samples' input size, plus (3 + number of samples) × the reference size when training, plus the model, plus 50 GB. | sized from inputs |

:::note One sample or many

- Pairing is by position. The run fails before any sample is processed if `gex_h5ads` is empty, `input_ids` differs in length, or `atac_h5ads` is neither empty nor the same length. ATAC is all or none per run.
- The label-transfer task then checks every input file (it exists, is not empty, opens as an h5ad) and that the ids are unique, and fails in seconds naming every problem, before any sample is labelled.
- In an inputs JSON a single sample still needs brackets: `"scANVI.input_ids": ["10k_PBMC"]`, `"scANVI.gex_h5ads": ["gs://…"]`. A Terra data-table column such as `this.gex_h5ad` maps to an array input directly, because Terra wraps the single value.
- Samples run one after another on one machine, so wall time is the sum of the samples. One bad sample fails the whole run, and a retry redoes every sample; call caching applies to the run as a whole.

:::

:::note GEX-only mode
The ATAC input is optional. Leave `atac_h5ads` empty (the default) and every sample runs in GEX-only mode: labels are transferred to GEX only, and `atac_annotated_h5ads` is empty. See [`example_inputs/scANVI.gex_only.json`](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/scanvi/example_inputs).
:::

#### Using a saved model

Every sample emits its SCANVI model in `scanvi_models_out`, a small `.tar.gz` with no bundled data. Pass a model as `scanvi_model` and every sample of the run is labelled with it, without training SCVI/SCANVI on the reference and without loading the reference at all. The pipeline reads the model's saved batch categories and chooses one of two paths per sample:

- **The model knows the sample's batch label.** A model this pipeline trained does, because the pipeline labels every sample's GEX batch `pd-multiome_sci_gex` (and ATAC `pd-multiome_sci_atac`). It is loaded and predicts as is.
- **The batch is new to the model.** A model trained on the reference alone knows only the reference's batches. The sample is mapped onto it with scArches: `SCANVI.load_query_data` followed by `query_max_epochs` of training on the sample alone, with the reference weights frozen.

Either way, the predictions file holds the sample's cells only. For a cheap run, also set `gpu_count = 0` and smaller `mem_size` / `nthreads`. Training once and then labelling many samples on CPU this way is far cheaper than training per sample (see [Example run costs](#example-run-costs)).

#### Running on Terra: per sample or per set

Terra can feed the sample arrays from a data table in two ways. The column names below (`sample_id`, `gex_h5ad`, `atac_h5ad`) are examples; use your own table's names.

- **Per row** (root entity `sample`): set `input_ids` = `this.sample_id`, `gex_h5ads` = `this.gex_h5ad`, `atac_h5ads` = `this.atac_h5ad`, and `scanvi_model` or `ref_h5ad` from the workspace. Each selected row runs as its own workflow on its own machine, and Terra wraps each single value into a one-element array. A row with no ATAC value runs GEX-only. Outputs land on the row as one-element lists.
- **Per set** (root entity `sample_set`): set `input_ids` = `this.samples.sample_id`, `gex_h5ads` = `this.samples.gex_h5ad` and `atac_h5ads` = `this.samples.atac_h5ad` (leave it unset for a GEX-only set). One workflow runs per set, on one machine. Outputs land on the set row in sample order. A set that mixes multiome and GEX-only samples fails the pairing check, so split it.

If a task runs out of memory, turn on Terra's **Retry with more memory** option when you launch the workflow: the task's one retry then runs with its memory multiplied.

Without a data table, list the samples in the inputs JSON. For two multiome samples labelled with a saved model (placeholder paths):

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
  "scANVI.scanvi_model": "gs://your-bucket/scanvi_model.tar.gz",
  "scANVI.gpu_count": 0
}
```

#### Reference requirements

When training, the reference h5ad must provide cell-type labels and a batch. Two reference forms are supported:

- **PBMC-style** (default): counts in `.X`, labels in `obs['final_annotation']`, batch in `obs['batch']`.
- **AIT-schema** (Allen Institute Taxonomy): auto-detected via `uns['schema_version']` + `uns['hierarchy']`. These typically have **no `.X`** (counts are materialized from `.raw`, keeping the gene symbols from `.var`), a cell-type **hierarchy** (e.g. `class`/`subclass`/`cluster_id`), and a batch column such as `donor_id`. By default the label is taken from `subclass` and the batch from `donor_id`; override with `ref_label_column` / `ref_batch_column`. For a mouse AIT reference used in multiome mode, set `genome` to `mm10`/`mm39` so the ATAC gene-activity conversion uses the correct genome.

The query datasets (GEX, and ATAC when present) do not need pre-existing annotations — placeholder `Unknown` labels are added automatically before training.

## scANVI tasks and tools

The [scANVI workflow](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.wdl) defines two tasks inline, plus `MultiomeLabelTransferCpu`, the CPU-only twin of MultiomeLabelTransfer. Both use the same Docker image, and only MultiomeLabelTransfer is allocated GPUs.

| Task name | Tool | Software | Description |
| --- | --- | --- | --- |
| [PreprocessReference](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.wdl) | scanpy | Python | Runs only when training, once per run: materializes counts from `.raw` for AIT references, maps the chosen label/batch columns to `final_annotation`/`batch`, and tags the reference cells with modality `rna_annotated`. |
| [MultiomeLabelTransfer](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.wdl) | scanpy, snapatac2, scvi-tools | Python (GPU, or CPU as MultiomeLabelTransferCpu) | Checks every input, then preprocesses and labels each sample in turn and writes its outputs. |

#### 1. PreprocessReference (CPU-only, training only)

Prepares the annotated reference once per run; a run with a supplied `scanvi_model` skips it. Steps:

1. **Load the reference** — Reads `ref_h5ad` (scanpy).
2. **Adapt the reference** — Detects an AIT-schema reference from `uns['schema_version']` plus `uns['hierarchy']`; materializes counts from `.raw` (keeping the `.var` gene symbols) when `.X` is absent; and maps the label and batch columns, from `ref_label_column` / `ref_batch_column` or the defaults (see [Reference requirements](#reference-requirements)), to `final_annotation` and `batch`.
3. **Modality tag** — Tags the reference cells `rna_annotated`.
4. **Reference id** — Names the reference (see [Model names](#model-names)); every model trained on it carries that name.
5. **Write output** — `preprocessed_ref.h5ad`, which the label-transfer task reads.

#### 2. MultiomeLabelTransfer (GPU or CPU)

Runs once per run. It first checks every input with `validate_inputs.py`, then runs the container's `scanvi_sample.py` once per sample, each in its own process and scratch directory, so memory and disk are freed before the next sample. For each sample:

1. **Preprocess** (`preprocess_query.py`) — Reads the sample's GEX (plus its ATAC cell-by-bin in multiome mode); adds `star_IsCell = True` and `gex_barcodes` when absent; keeps STARsolo cell calls; removes genes and cells with fewer than 3 counts. In multiome mode it reindexes ATAC by `gex_barcodes`, keeps the barcodes shared with GEX, and converts the cell-by-bin matrix into a gene activity matrix with `snapatac2.pp.make_gene_matrix` for `genome`. It sets the batch labels (`pd-multiome_sci_gex`, `pd-multiome_sci_atac`), the placeholder label `Unknown`, and the modality tags.
2. **Obtain the model** — Trains SCVI then SCANVI on the reference plus the sample when no model is supplied (`run_multi_model` / `run_gex_only_model`: genes in ≥ 5 cells, 5,000 highly variable genes, a 2-layer SCVI with 30 latent dimensions, up to 500 epochs each). With a supplied model, loads it or maps the sample onto it (see [Using a saved model](#using-a-saved-model)).
3. **Predict labels** — `transfer_labels()` predicts a cell type (`C_scANVI`) for every cell and extracts the scANVI latent (`X_scANVI`).
4. **Write the annotated matrices** — Copies each cell's label (and `max_probability` when asked) and latent row into the GEX (and ATAC) matrix, then makes it SCP-ready (`scp_output.py`; see [Outputs](#outputs)). Writes `~{input_id}_gex_annotated_matrix.h5ad` (and `~{input_id}_atac_annotated_matrix.h5ad`).
5. **Write the predictions and the model** — `finalize_output()` renames `final_annotation` → `celltype` and copies the counts into `.raw`; writes `~{input_id}_SCANVI_predictions.h5ad` and `~{input_id}_scanvi_model.tar.gz`.

The task lists each sample's output files in input order as it goes and reads its output arrays from those lists, so `outputs[i]` always belongs to `input_ids[i]`.

### Workflow diagram

```text
  ref_h5ad (training only)          input_ids, gex_h5ads (+ atac_h5ads), scanvi_model
         │                                       │
         ▼                                       │
┌───────────────────────┐                        │
│  PreprocessReference  │  CPU, once per run     │
└───────────┬───────────┘                        │
            │ preprocessed_ref.h5ad              │
            ▼                                    ▼
┌───────────────────────────────────────────────────────────────┐
│                   MultiomeLabelTransfer                       │  one machine
│        (MultiomeLabelTransferCpu when gpu_count = 0)          │
│                                                               │
│  check every input                                            │
│  for each sample i, one after another:                        │
│      preprocess → obtain model → predict → write outputs      │
└───────────────────────────────┬───────────────────────────────┘
                                │
                                ▼
           input_ids[i]_gex_annotated_matrix.h5ad    (SCP-ready)
           input_ids[i]_atac_annotated_matrix.h5ad   (multiome only)
           input_ids[i]_SCANVI_predictions.h5ad
           input_ids[i]_scanvi_model.tar.gz
```

### Design rationale

- **One machine per run.** On Terra, starting a machine per sample costs more in waiting than in work for typical samples: the image is 11 GB, and large machine shapes hit zone-capacity retries. Running samples one after another on one machine starts it once; each sample's own process frees its memory, so the machine only has to fit the largest sample.
- **No reference with a supplied model.** Predicting, embedding and laying out the sample alone keeps a run's cost independent of the reference's size, which matters for whole-brain references of millions of cells.
- **GPU only where it helps.** The GPU and CPU twins run the same command; scvi-tools picks the accelerator itself, so a supplied-model run can use the cheaper CPU machine.

### Outputs

Every output except `pipeline_version_out` is an array with one entry per sample, in the order of `input_ids`; each file name starts with that sample's input id.

| Output Variable Name | Filename | Output Type | Output Format |
| --- | --- | --- | --- |
| `gex_annotated_h5ads` | `<input_id>_gex_annotated_matrix.h5ad` | The sample's filtered GEX cells with their transferred labels (`final_annotation`), SCP-ready (below). | H5AD |
| `atac_annotated_h5ads` | `<input_id>_atac_annotated_matrix.h5ad` | The sample's ATAC gene-activity cells with their transferred labels, in the same form. Produced only in multiome mode; empty in GEX-only mode. | H5AD |
| `scanvi_predictions_h5ads` | `<input_id>_SCANVI_predictions.h5ad` | The labelled data the model saw, with SCANVI predictions (`celltype`), latent and UMAP; counts in `.raw`. The sample's cells only with a supplied model; the sample plus the reference when training. Its `obs` keeps only the columns shared by every dataset combined into it (GEX, ATAC, and the reference when training), so a column only one input has, such as a reference-only `donor_id`, is dropped; use the annotated matrices for each input's own columns. | H5AD |
| `scanvi_models_out` | `<input_id>_scanvi_model.tar.gz` | The sample's trained, loaded or mapped SCANVI model, saved without bundled AnnData, with its reference's name in `reference.txt`. Can be passed back as `scanvi_model`. | TAR.GZ |
| `pipeline_version_out` | N/A | Version of the processing pipeline run on this data. | String |

The annotated matrices are written in the form SCP reads, compressed with gzip:

| Slot | Contents |
| --- | --- |
| `.raw` | The counts (float32), over every gene the matrix keeps. |
| `.X` | log2(CPM + 1) (float32): counts per million, then log base 2, so a difference of 1 is a doubling and the units are bits. This is the Allen Institute's convention, the same as the ABC Atlas log2 matrices and the `log2CPM` input of MapMyCells; `uns['log1p']['base']` records the base as 2. |
| `obsm['X_umap']` | A 2-column UMAP fitted on the sample's own latent: the layout SCP plots. Each sample fits its own, so UMAP axes are not comparable between samples. |
| `obsm['X_scANVI']` | The full scANVI latent (30 dimensions). Comparable between samples labelled with the same model. |
| `obs['final_annotation']` | The transferred cell-type label (plus `max_probability` when requested). |
| `uns` | Run provenance: `scanvi_input_id`, `scanvi_input_file` (the GEX file's name) and `scanvi_model`, the model's name: the reference it was trained on (see [Model names](#model-names)). |

#### Model names

A model is named by the reference it was trained on: the reference's `uns['title']` (the Allen taxonomy and CELLxGENE schemas both define it), or its file name when it has none, made file-name safe, followed by the first 8 hex digits of the reference file's sha256. For example, `Mouse_cortex_hippocampus_10X_seq_04042025-f3177e00`. The title says which reference; the hash tells apart subsets and rebuilds that share a title.

A training run computes the name in PreprocessReference and writes it into each model archive as `reference.txt`, next to `model.pt`; a supplied model's `reference.txt` is read and carried forward. Every output records it as `uns['scanvi_model']`. A model archive without `reference.txt`, such as one from before scANVI 3.0.0, is named by its archive's file name instead. Each model archive holds one SCANVI model; the scVI encoder and decoder it was built from are part of it, but there is no separate scVI model.

:::note Run metadata, not sample metadata
The annotated matrices keep every `obs` column of their own input (the Optimus QC columns, for example) plus the transferred label; the predictions file keeps only columns shared by the datasets combined into it. Each output records which sample, input file and model produced it, in `uns` under keys prefixed `scanvi_`; the prefix keeps them apart from the `uns['input_id']` and `uns['pipeline_version']` an Optimus h5ad already carries. The pipeline adds no sample-level metadata (donor, species, organ, disease, library preparation, sex) to any output. Sample metadata is joined downstream from the data archive's sample manifests, not by this pipeline.
:::

## Runtime configuration

All tasks use the same Docker image (pinned by digest). GPU and CUDA setup is handled entirely by the execution engine — the container does not configure the GPU environment itself.

#### Task 1 — `PreprocessReference` (CPU-only, once per run, training only)

| Attribute | Value |
| --- | --- |
| `docker` | `us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:5764b99116daedef6c73a83483028ac78e318b1ec33de8824a9861b4d3201b52` |
| `bootDiskSizeGb` | 20 |
| `disks` | `local-disk` 3 × the reference size + 20 GB, SSD |
| `memory` | `104 GiB` |
| `cpu` | 16 |
| `maxRetries` | 1 |

#### Task 2 — `MultiomeLabelTransfer` (GPU; CPU-only as `MultiomeLabelTransferCpu`)

| Attribute | Value |
| --- | --- |
| `docker` | `us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:5764b99116daedef6c73a83483028ac78e318b1ec33de8824a9861b4d3201b52` |
| `bootDiskSizeGb` | 20 |
| `disks` | `local-disk` `disk_size`, or sized from the inputs (see [Workflow inputs](#workflow-inputs)), SSD |
| `memory` | `mem_size` (default `104 GiB`) |
| `cpu` | `nthreads` (default 16) |
| `gpuType` | `nvidia-tesla-t4` (GPU task only) |
| `gpuCount` | `gpu_count` (default 2; GPU task only) |
| `nvidiaDriverVersion` | `535.104.05` (GPU task only) |
| `maxRetries` | 1 |

:::note GPU driver compatibility
Driver version `535.104.05` is compatible with CUDA 12.x and NVIDIA T4 GPUs and has been verified working on GCP / Terra with the `scvi-scanvi` container.
:::

:::note CPU-only variant
When `gpu_count = 0` the workflow routes to `MultiomeLabelTransferCpu` — an identical task with **no** GPU runtime attributes (Cromwell rejects `gpuCount = 0` and can't conditionally omit the GPU keys from one task). Both tasks run the same command; scvi-tools auto-detects the accelerator, so it uses the GPU when present and CPU otherwise.
:::

## Example run costs

The table below reports representative Terra costs for scANVI runs across three references — mouse hippocampus, human neocortex, and PBMC — in both training and CPU-inference modes. Costs are from single runs on GCP `us-central1` (training on **2× NVIDIA T4** GPUs; inference CPU-only via `gpu_count = 0`) and will vary with region, machine type, preemption, and dataset. "Cells labelled" is the number of query cells annotated (after preprocessing; for multiome, after the GEX↔ATAC shared-barcode intersection); "Reference" is the annotated `ref_h5ad`.

Each row is a single-sample run measured before the one-machine layout, when each sample ran on its own machines. A run of many samples now starts one machine and labels the samples in turn. Without `scanvi_model` each sample also trains its own model, so a many-sample training run costs roughly the train cost for that reference times the number of samples. Labelling every sample with one saved model on CPU (`scanvi_model` plus `gpu_count = 0`) is far cheaper.

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
