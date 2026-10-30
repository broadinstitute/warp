---
sidebar_position: 1
slug: /Pipelines/scANVI_Pipeline/README
---

# scANVI Overview

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| See [changelog](https://github.com/broadinstitute/warp/blob/develop/pipelines/wdl/scanvi/scANVI.changelog.md) for version information. | See changelog | WARP Pipelines | Please [file an issue in WARP](https://github.com/broadinstitute/warp/issues)  |

## Introduction to the scANVI workflow

The scANVI pipeline is a cloud-optimized WDL workflow that performs **cell type label transfer** using [scVI](https://docs.scvi-tools.org/en/stable/user_guide/models/scvi.html) and [scANVI](https://docs.scvi-tools.org/en/stable/user_guide/models/scanvi.html) (single-cell ANnotation using Variational Inference) deep generative models. It labels single-cell RNA-seq (GEX) data and, optionally, the matching ATAC-seq data (multiome mode) with the cell types of an annotated reference.

One run labels **one or more samples**, each independently and with its own outputs, on **one machine, one after another**. Each sample is labelled at one or more **label levels** of the reference, such as `class` and `subclass`, with one SCANVI model per level. Each input file `<name>.h5ad` gives one `<name>_annotated.h5ad`, ready for [Single Cell Portal](https://singlecell.broadinstitute.org) (SCP), which records which sample, input file, model and levels produced it.

The models are either trained in the run, on the reference plus the sample, or supplied as a saved model archive; a sample is mapped onto a model trained on the reference alone with scArches (see [How the label transfer works](#how-the-label-transfer-works)).

:::tip Want to use scANVI for your publication?
The pipeline is designed to consume the outputs of the [Multiome](../Multiome_Pipeline/README.md) and [PeakCalling](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/peak_calling) WARP pipelines. Identify it in your methods section with its SciCrunch resource identifier, *scANVI Pipeline (RRID:SCR_028705)*, and cite the pipeline using the WARP citation in the [Citing](#citing-the-scanvi-pipeline) section below.
:::

### How the label transfer works

scANVI annotates in two stages, an **unsupervised** representation-learning step followed by a **semi-supervised** annotation step:

1. **SCVI, unsupervised embedding.** SCVI is a variational autoencoder trained **without cell-type labels**. It learns a low-dimensional latent space that captures biological variation while correcting batch effects across donors, between the sample and the reference, and across modalities.
2. **SCANVI, semi-supervised annotation.** SCANVI is initialized **from the trained SCVI model** and adds a cell-type classifier on that latent space. Reference cells supply their known labels; unlabelled cells (tagged `Unknown`) keep shaping the latent space. SCANVI then predicts a cell type for every cell. Each label level gets its own SCANVI model from the same SCVI model, so a run with the levels `class` and `subclass` predicts both for every cell.

The models for a run come from one of three places, chosen per run (and per level for a supplied model):

| Models | What happens to each sample | Reference needed |
| --- | --- | --- |
| None supplied (`ref_h5ad` given) | SCVI is trained on the reference plus the sample, then a SCANVI per level; the models are saved in `scanvi_models_out`. | Yes |
| A model this pipeline trained | The model already knows the sample's batch label, so each level is loaded and predicts as is. | No |
| A model trained on the reference alone, e.g. by [warp-scanvi-training](https://github.com/rsc3/warp-scanvi-training) | The sample's batch is new to the model, so the sample is mapped onto each level with scArches (below). | No |

A supplied model labels every sample of the run, so training once and labelling many samples, on CPU if need be, is far cheaper than training per sample (see [Example run costs](#example-run-costs)).

### Query mapping with scArches

scArches ([Lotfollahi et al. 2022](https://doi.org/10.1038/s41587-021-01001-7)) maps new data onto a trained model without retraining it. The pipeline uses it when a supplied model does not know the sample's batch label, which is the case for every model trained on a reference alone.

**What the pipeline does.** For each level, `SCANVI.load_query_data` adds the sample's batch label to the model, and the model trains for `query_max_epochs` (default 100, the scvi-tools tutorial value) on the sample alone. scvi-tools' defaults keep the reference's weights fixed, the classifier included, and train only the weights for the new batch in the first layers of the encoder and decoder; weight decay is off, so the fixed weights stay exactly as trained. The sample's cells then sit in the reference's latent space, and the reference's classifier labels them. Each sample is mapped in its own process, on its own cells; the reference is never loaded.

**What the model needs.** Mapping works best on a model trained with scvi-tools' scArches settings: the batch covariate encoded (`encode_covariates=True`), layer norm instead of batch norm, and gene-level dispersion, because a new batch's per-batch dispersion would stay untrained. warp-scanvi-training uses these settings. Models trained by this pipeline do not, but they are never mapped: they already know the pipeline's batch label.

**What comes out.** Each level gives the sample's labels, its latent in the reference's space (comparable between samples mapped onto the same model) and a UMAP of the sample's own cells. Mapping time grows with the number of levels, since each level is mapped separately.

**Status and open questions.**

- **Tested on the Plumbing data only.** Mapped as a new sample, a reference's own cells get their known labels back (98.8 % on the 800-cell Plumbing reference). Accuracy on real samples has not been measured yet (coming soon).
- **Run it on a GPU for large references.** A warp-scanvi-training model covers every gene of its atlas, and mapping onto such a model is slow on CPU (about 8 s per epoch for 800 cells on 4 CPUs). Set `gpu_count` to 1 or more for real samples.
- **Levels are mapped independently**, so a cell's labels at two levels can disagree with the reference's hierarchy: about 9 % of cells for class versus subclass on the Plumbing data, with tiny models and 20 mapping epochs. How often this happens on real data is not yet measured.
- **Multiome samples on a reference-only model run but are untested.** The ATAC gene-activity matrix is a stand-in for expression that scVI's assumptions fit poorly.
- **One batch label per sample.** Every sample's GEX is labelled `pd-multiome_sci_gex` (ATAC `pd-multiome_sci_atac`). On the mapping path this is harmless, because each sample is mapped on its own; on a pipeline-trained model it means every sample shares the training sample's batch, so sample-to-sample batch effects are not corrected there.

:::note Coming soon

- **Validated query mapping.** Accuracy on real samples, starting with the Scientific mouse hippocampus data and a warp-scanvi-training model.
- **The mouse whole brain model.** The AIT33 (BICCN whole mouse brain consensus) model is being retrained with warp-scanvi-training.

:::

### ATAC gene-activity conversion (multiome mode)

SCVI needs every dataset it combines to share one feature space, so the ATAC cell-by-bin matrix cannot be combined with the cell-by-gene RNA matrices directly. The pipeline projects the ATAC signal onto genes instead:

1. **No peak calling.** Rather than relying on distal peaks, which are hard to assign to genes, the pipeline counts the Tn5 insertion fragments that overlap each gene's regulatory domain.
2. **Gene activity scoring.** Using the genome annotation for `genome` (hg38, mm10 or mm39), `snapatac2.pp.make_gene_matrix` takes each gene's body plus a promoter window upstream of its transcription start site and counts the fragments in it; fragments in distant intergenic regions are ignored.
3. **Combination.** The result is a cell-by-gene "activity" matrix that resembles an RNA count matrix, so the GEX, the ATAC activity and (when training) the reference can share one SCVI latent space, with the modality handled as a batch effect.

### How the models are trained

Training and mapping use **minibatch stochastic gradient descent**: the data are held in host memory, and each step sends a minibatch of cells (`batch_size`, default 128) to the GPU, updates the weights and fetches the next. GPU memory therefore scales with *minibatch size × number of genes*, not with the number of cells, so a large sample can be labelled on a modest GPU while the data live in host RAM.

**Two meanings of "batch."** `batch_size` is that SGD minibatch, an optimization detail. The **batch covariate** is the experimental grouping (donor, library, modality) that SCVI/SCANVI model in order to remove it as technical variation. The two are unrelated.

## Quickstart table

| Pipeline features | Description | Source |
| --- | --- | --- |
| Assay type | 10x single-cell / single-nucleus GEX, or Multiome (GEX + ATAC) | [10x Genomics](https://www.10xgenomics.com) |
| Overall workflow | Reference preparation when training (CPU, once per run), then every sample preprocessed and labelled in turn on one GPU or CPU machine | Code available on [GitHub](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.wdl) |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Models | SCVI (unsupervised VAE) + one SCANVI (semi-supervised classifier) per label level; scArches query mapping | [scvi-tools 1.5.1](https://docs.scvi-tools.org/) |
| ATAC gene-activity conversion | Cell-by-bin matrix → gene activity matrix (hg38, mm10 or mm39) | [snapatac2 2.10](https://kzhang.org/SnapATAC2/) |
| Data input format | AnnData h5ad files: one or more GEX samples, optionally one ATAC cell-by-bin file per sample; plus a model archive or an annotated reference | [AnnData](https://anndata.readthedocs.io/) |
| Data output format | One annotated h5ad per input file, with labels, latent and UMAP per label level | [AnnData](https://anndata.readthedocs.io/) |

## Set-up

### scANVI installation

To download the latest scANVI release, see the release tags prefixed with "scANVI" on the WARP [releases page](https://github.com/broadinstitute/warp/releases). All scANVI pipeline releases are documented in the [scANVI changelog](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.changelog.md).

scANVI can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant, flexible workflow management system that supports multiple computing platforms. The workflow can also be run in [Terra](https://app.terra.bio).

### Inputs

scANVI takes its samples as arrays paired by position: `input_ids[i]`, `gex_h5ads[i]` and, in multiome mode, `atac_h5ads[i]` describe sample *i*. Outputs are named after the input files, so the input files of a run must have distinct names. Give either `scanvi_model` (to label with a saved model) or `ref_h5ad` (to train); the run fails at once with neither.

Example input JSON files, including a GEX-only one, are in the [`example_inputs`](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/scanvi/example_inputs) folder.

#### Workflow inputs

| Input | Type | Description | Default |
| --- | --- | --- | --- |
| `input_ids` | Array[String] | One identifier per sample (for example the library's id), in the same order as `gex_h5ads`. Recorded in each of the sample's outputs as `uns['scanvi_input_id']`; ids must be unique. | — (required) |
| `gex_h5ads` | Array[File] | Gene expression AnnData h5ad files from Multiome / Optimus, one per sample. | — (required) |
| `atac_h5ads` | Array[File] | **Optional** ATAC cell-by-bin AnnData h5ad files from Multiome / PeakCalling, one per `gex_h5ads` entry and in the same order. Empty runs every sample in GEX-only mode, and `atac_annotated_h5ads` is then empty. | `[]` |
| `scanvi_model` | File? | **Optional** model archive (see [Models](#models)), used for every sample instead of training; every level labels every sample. | — |
| `ref_h5ad` | File? | Annotated reference h5ad (see [Reference requirements](#reference-requirements)). Needed only to train, that is, when `scanvi_model` is not given; ignored otherwise. | — |
| `ref_label_columns` | Array[String] | When training, the label levels: the reference `obs` columns to label each sample with, one SCANVI model each, for example `["class", "subclass"]`. The first gives the default UMAP. When empty, defaults to `["subclass"]` for AIT references and `["final_annotation"]` otherwise. | `[]` |
| `ref_batch_column` | String? | When training, the reference `obs` column to use as the batch. When unset, defaults to `donor_id` for AIT references and `batch` otherwise. | — |
| `max_epochs` | Int? | Cap on SCVI/SCANVI training epochs when training. When unset, the container default (500) is used. | — |
| `query_max_epochs` | Int? | Epochs when mapping a sample onto a reference-only model (see [Query mapping with scArches](#query-mapping-with-scarches)). When unset, the container default (100) is used. | — |
| `batch_size` | Int | SCVI/SCANVI minibatch size. Lower it to fit a high-cardinality reference on a small GPU (activation memory scales with `batch_size` × number of labels); raise it on a large-VRAM GPU. | `128` |
| `seed` | Int? | Random seed for SCVI/SCANVI training and mapping. Set it to make a run reproducible (the test wrappers pin one); unset keeps training nondeterministic. | — |
| `genome` | String | Genome for the ATAC gene-activity conversion (multiome only): `hg38`, `mm10` or `mm39`. | `"hg38"` |
| `output_max_probability` | Boolean | Add a `<level>_max_probability` obs column per level (the per-cell maximum SCANVI posterior: the assigned label's confidence) to every annotated h5ad. Set `false` to leave it out. | `true` |
| `gpu_count` | Int | GPUs for the label-transfer machine. `0` gives a CPU-only run, which suits a supplied model that needs no mapping. | `2` |
| `mem_size` | Int | Memory (GiB) for the label-transfer machine. Samples run one at a time, so it only has to fit the largest sample (plus the reference when training). | `104` |
| `nthreads` | Int | CPUs for the label-transfer machine. | `16` |
| `disk_size` | Int? | Disk (GB) for the label-transfer machine. When unset: 4 × the samples' input size, plus 3 × the reference size when training, plus the model, plus 50 GB. | sized from inputs |

:::note One sample or many

- Pairing is by position. The run fails before any sample is processed if `gex_h5ads` is empty, `input_ids` differs in length, or `atac_h5ads` is neither empty nor the same length. ATAC is all or none per run.
- The label-transfer task then checks every input file (it exists, is not empty, opens as an h5ad), that the ids are unique and that no two input files share a name, and fails in seconds naming every problem, before any sample is labelled.
- In an inputs JSON a single sample still needs brackets: `"scANVI.input_ids": ["10k_PBMC"]`, `"scANVI.gex_h5ads": ["gs://…"]`. A Terra data-table column such as `this.gex_h5ad` maps to an array input directly, because Terra wraps the single value.
- Samples run one after another on one machine, so wall time is the sum of the samples. One bad sample fails the whole run, and a retry redoes every sample; call caching applies to the run as a whole.

:::

#### Running on Terra: per sample or per set

Terra can feed the sample arrays from a data table in two ways. The column names below (`sample_id`, `gex_h5ad`, `atac_h5ad`) are examples; use your own table's names.

- **Per row** (root entity `sample`): set `input_ids` = `this.sample_id`, `gex_h5ads` = `this.gex_h5ad`, `atac_h5ads` = `this.atac_h5ad`, and `scanvi_model` or `ref_h5ad` from the workspace. Each selected row runs as its own workflow on its own machine, and Terra wraps each single value into a one-element array. A row with no ATAC value runs GEX-only. Outputs land on the row as one-element lists.
- **Per set** (root entity `sample_set`): set `input_ids` = `this.samples.sample_id`, `gex_h5ads` = `this.samples.gex_h5ad` and `atac_h5ads` = `this.samples.atac_h5ad` (leave it unset for a GEX-only set). One workflow runs per set, on one machine. Outputs land on the set row in sample order. A set that mixes multiome and GEX-only samples fails the pairing check, so split it.

If a task runs out of memory, turn on Terra's **Retry with more memory** option when you launch the workflow: the task's one retry then runs with its memory multiplied.

Without a data table, list the samples in the inputs JSON. For two GEX samples mapped onto a saved model on one GPU (placeholder paths):

```json
{
  "scANVI.input_ids": ["sample_a", "sample_b"],
  "scANVI.gex_h5ads": [
    "gs://your-bucket/sample_a_gex.h5ad",
    "gs://your-bucket/sample_b_gex.h5ad"
  ],
  "scANVI.scanvi_model": "gs://your-bucket/reference_scanvi_model.tar.gz",
  "scANVI.gpu_count": 1
}
```

#### Reference requirements

When training, the reference h5ad must provide cell-type labels and a batch. Two reference forms are supported:

- **PBMC-style** (default): counts in `.X`, labels in `obs['final_annotation']`, batch in `obs['batch']`.
- **AIT-schema** (Allen Institute Taxonomy): auto-detected via `uns['schema_version']` + `uns['hierarchy']`. These typically have **no `.X`** (counts are materialized from `.raw`, keeping the gene symbols from `.var`), a cell-type hierarchy (e.g. `class`/`subclass`/`cluster_id`), and a batch column such as `donor_id`. By default the label is taken from `subclass` and the batch from `donor_id`; override with `ref_label_columns` (one or more levels) and `ref_batch_column`. For a mouse reference in multiome mode, set `genome` to `mm10` or `mm39`.

The samples need no annotations; they are labelled `Unknown` (unlabelled) for training and mapping.

## Models

### Model archives

A model archive is one `.tar.gz` holding every model trained from one reference, so a reference's models are tracked as a single file:

```text
scanvi_model/
  reference.txt             the model's name (below)
  levels.txt                the label levels, one per line, in order; the first gives the default UMAP
  scvi/model.pt             the scVI model the SCANVI models were trained from
  scanvi/<level>/model.pt   one SCANVI model per label level
```

Every model is saved without data. Level names become `obs` columns and part of `obsm` keys, so they are limited to letters, digits, `_`, `.` and `-`. The archive is written and read by the image's `model_archive.py`; an archive in the older single-model form is refused with a message naming this layout.

A sample's genes are matched to the model's by name. A model whose genes are Ensembl IDs, such as a warp-scanvi-training model trained with `--gene-names ensembl_id`, is matched through the sample's Optimus `ensembl_ids` column instead, which finds genes renamed between annotation versions: a SCORCH sample has 99.5 % of the AIT33 atlas's genes by Ensembl ID and 93.3 % by symbol. The annotated outputs keep the sample's own gene names, and the log reports how many of the model's genes the sample has.

Archives come from two places:

- **A training run of this pipeline:** `scanvi_models_out` holds one archive per sample, trained on the reference plus that sample. Passed back as `scanvi_model`, it loads as is.
- **[warp-scanvi-training](https://github.com/rsc3/warp-scanvi-training):** trains scVI and a SCANVI per level on a whole atlas out of core, with the scArches settings, and writes the same archive with the same `model_archive.py`. Samples are mapped onto it (see [Query mapping with scArches](#query-mapping-with-scarches)).

### Model names

A model is named by the reference it was trained on: the reference's `uns['title']` (the Allen taxonomy and CELLxGENE schemas both define it), or its file name when it has none, made file-name safe, followed by the first 8 hex digits of the reference file's sha256. For example, `Mouse_cortex_hippocampus_10X_seq_04042025-f3177e00`. The title says which reference; the hash tells apart subsets and rebuilds that share a title.

A training run computes the name in PreprocessReference and writes it into the archive as `reference.txt`; warp-scanvi-training does the same by the same rule. Every output records the name as `uns['scanvi_model']`.

## scANVI tasks and tools

The [scANVI workflow](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.wdl) defines two tasks inline, plus `MultiomeLabelTransferCpu`, the CPU-only twin of MultiomeLabelTransfer. Both use the same Docker image, and only MultiomeLabelTransfer is allocated GPUs.

| Task name | Tool | Software | Runs |
| --- | --- | --- | --- |
| [PreprocessReference](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.wdl) | scanpy | Python, CPU | Once per training run |
| [MultiomeLabelTransfer](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.wdl) | scanpy, snapatac2, scvi-tools | Python, GPU (or CPU as MultiomeLabelTransferCpu) | Once per run, every sample in turn |

### 1. PreprocessReference (CPU-only, training only)

Prepares the annotated reference once per run; a run with a supplied `scanvi_model` skips it.

1. **Adapt the reference** — Reads `ref_h5ad`; detects an AIT-schema reference from `uns['schema_version']` plus `uns['hierarchy']`; materializes counts from `.raw` (keeping the `.var` gene symbols) when `.X` is absent; keeps the label levels from `ref_label_columns` (or the defaults) as string columns, listed in `uns['scanvi_levels']`; and maps the batch column to `batch`.
2. **Reference id** — Names the reference (see [Model names](#model-names)).
3. **Write output** — `preprocessed_ref.h5ad`, which the label-transfer task reads.

### 2. MultiomeLabelTransfer (GPU or CPU)

Checks every input with `validate_inputs.py`, then runs the container's `scanvi_sample.py` once per sample, each in its own process and scratch directory, so memory and disk are freed before the next sample. For each sample:

1. **Preprocess** (`preprocess_query.py`) — Reads the sample's GEX (plus its ATAC cell-by-bin in multiome mode); adds `star_IsCell = True` and `gex_barcodes` when absent; keeps STARsolo cell calls; removes genes and cells with fewer than 3 counts. In multiome mode it reindexes ATAC by `gex_barcodes`, keeps the barcodes shared with GEX, and converts the cell-by-bin matrix into a gene activity matrix. It sets the batch labels (`pd-multiome_sci_gex`, `pd-multiome_sci_atac`).
2. **Obtain the models and predict** — When training, trains SCVI once on the reference plus the sample (genes in ≥ 5 cells, 5,000 highly variable genes, 2 layers, 30 latent dimensions, up to `max_epochs`), then a SCANVI from it per level, and saves them all as one archive. With a supplied model, loads each level or maps the sample onto it. Each level's model predicts a label for every cell and gives its latent.
3. **Write the annotated matrices** — Copies each cell's labels and latent rows into the GEX (and ATAC) matrix, makes it SCP-ready (`scp_output.py`) and writes it (see [Outputs](#outputs)).

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
│      preprocess → for each level: obtain model → predict      │
│      → write outputs                                          │
└───────────────────────────────┬───────────────────────────────┘
                                │
                                ▼
           <gex file name>_annotated.h5ad       (SCP-ready)
           <atac file name>_annotated.h5ad      (multiome only)
           <gex file name>_scanvi_model.tar.gz  (training only)
```

### Design rationale

- **One machine per run.** On Terra, starting a machine per sample costs more in waiting than in work for typical samples: the image is 11 GB, and large machine shapes hit zone-capacity retries. Each sample's own process frees its memory, so the machine only has to fit the largest sample.
- **No reference with a supplied model.** Labelling and laying out the sample alone keeps a run's cost independent of the reference's size, which matters for whole-brain references of millions of cells.
- **One archive per reference.** A reference's scVI and SCANVI models are small and belong together, so they travel as one file, named by the reference.
- **GPU only where it helps.** The GPU and CPU twins run the same command; scvi-tools picks the accelerator itself.

## Outputs

Every output except `pipeline_version_out` is an array with one entry per sample, in the order of `input_ids`. The workflow builds the file names from the input files before the task runs, so `outputs[i]` always belongs to `input_ids[i]`, and a missing output fails the task.

| Output Variable Name | Filename | Output Type | Output Format |
| --- | --- | --- | --- |
| `gex_annotated_h5ads` | `<gex file name>_annotated.h5ad` | The sample's filtered GEX cells with their labels at each level, SCP-ready (below). | H5AD |
| `atac_annotated_h5ads` | `<atac file name>_annotated.h5ad` | The sample's ATAC gene-activity cells with their labels, in the same form. Multiome mode only; empty in GEX-only mode. | H5AD |
| `scanvi_models_out` | `<gex file name>_scanvi_model.tar.gz` | Training runs only (empty with a supplied model): the sample's model archive (see [Models](#models)). | TAR.GZ |
| `pipeline_version_out` | N/A | Version of the processing pipeline run on this data. | String |

The annotated matrices are written in the form SCP reads, compressed with gzip:

| Slot | Contents |
| --- | --- |
| `.raw` | The counts (float32), over every gene the matrix keeps. |
| `.X` | log2(CPM + 1) (float32): counts per million, then log base 2, so a difference of 1 is a doubling and the units are bits. This is the Allen Institute's convention, the same as the ABC Atlas log2 matrices and the `log2CPM` input of MapMyCells; `uns['log1p']['base']` records the base as 2. |
| `obs['<level>']` | The label at each level, in a column named after the level, for example `obs['class']` and `obs['subclass']`, and next to each its confidence, `obs['<level>_max_probability']` (the maximum SCANVI posterior; left out when `output_max_probability` is `false`). A column of the input with the same name is replaced. |
| `obsm['X_scVI_<level>']` | Each level's 30-dimensional scVI latent: the encoder's latent variable z, as trained with that level's SCANVI model. Comparable between samples labelled with the same model. |
| `obsm['X_umap_<level>']` | A 2-column UMAP fitted on that level's latent, for the sample's own cells. Each sample fits its own, so UMAP axes are not comparable between samples. |
| `obsm['X_umap']` | A copy of the first level's UMAP: the layout SCP plots by default. |
| `uns` | Run provenance: `scanvi_input_id`, `scanvi_input_file` (the GEX file's name), `scanvi_model` (the model's name, see [Model names](#model-names)), `scanvi_levels` (the levels in order) and `scanvi_versions` (the key library versions of the image that wrote the file: scvi-tools, torch, scanpy, anndata, snapatac2, numpy, scipy). |

:::note Run metadata, not sample metadata
The annotated matrices keep every `obs` column of their own input (the Optimus QC columns, for example) plus the labels and the few columns preprocessing adds (see [Annotated h5ad spec](#annotated-h5ad-spec)). The `scanvi_` prefix on the provenance keys keeps them apart from the `uns['input_id']` and `uns['pipeline_version']` an Optimus h5ad already carries. The pipeline adds no sample-level metadata (donor, species, organ, disease, library preparation, sex); sample metadata is joined downstream from the data archive's sample manifests.
:::

## Runtime configuration

All tasks use the same Docker image, pinned by digest. GPU and CUDA setup is handled by the execution engine; the container does not configure the GPU environment itself.

### Task 1 — `PreprocessReference` (CPU-only, once per run, training only)

| Attribute | Value |
| --- | --- |
| `docker` | `us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:938baa80261c9d2122838403d36e000cf914035705f4cbddc6eed34595289386` |
| `bootDiskSizeGb` | 20 |
| `disks` | `local-disk` 3 × the reference size + 20 GB, SSD |
| `memory` | `104 GiB` |
| `cpu` | 16 |
| `maxRetries` | 1 |

### Task 2 — `MultiomeLabelTransfer` (GPU; CPU-only as `MultiomeLabelTransferCpu`)

| Attribute | Value |
| --- | --- |
| `docker` | `us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:938baa80261c9d2122838403d36e000cf914035705f4cbddc6eed34595289386` |
| `bootDiskSizeGb` | 20 |
| `disks` | `local-disk` `disk_size`, or sized from the inputs (see [Workflow inputs](#workflow-inputs)), SSD |
| `memory` | `mem_size` (default `104 GiB`) |
| `cpu` | `nthreads` (default 16) |
| `gpuType` | `nvidia-tesla-t4` (GPU task only) |
| `gpuCount` | `gpu_count` (default 2; GPU task only) |
| `nvidiaDriverVersion` | `535.104.05` (GPU task only; compatible with CUDA 12.x and T4 GPUs, verified on GCP / Terra) |
| `maxRetries` | 1 |

:::note CPU-only variant
When `gpu_count = 0` the workflow routes to `MultiomeLabelTransferCpu`, an identical task with **no** GPU runtime attributes: Cromwell rejects `gpuCount = 0` and cannot conditionally omit the GPU keys from one task.
:::

## Example run costs

The table below reports representative Terra costs across three references (mouse hippocampus, human neocortex and PBMC), in training and CPU-inference modes, from single runs on GCP `us-central1` (training on **2× NVIDIA T4** GPUs; inference with `gpu_count = 0`). "Cells labelled" counts the sample's cells after preprocessing (for multiome, after the GEX↔ATAC shared-barcode intersection); "Reference" is the annotated `ref_h5ad`.

These runs predate the one-machine layout, label levels and scArches mapping: each labelled one sample at one level, and inference loaded a pipeline-trained model without mapping. A run of many samples now starts one machine; each extra level adds one SCANVI training or mapping per sample; mapping onto a large reference-only model needs a GPU. New measurements are planned.

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

- **Inference is far cheaper than training.** Loading a saved model and predicting on CPU cost 8–47× less than the equivalent training run.
- **Training cost tracks the reference, not the query.** Training cost scales with reference size (mouse 250,734 cells / 12.87 GB → $26.34 vs human 47,432 / 4.97 GB → $2.51).
- **Per-cell cost rises as the query shrinks**, because the fixed reference-load and training cost is spread over fewer cells (PBMC's ~900-cell queries cost the most per cell).
- **Multiome and GEX-only cost about the same** ($3.64 vs $3.38).

## Docker image

The `scvi-scanvi` image is maintained in [warp-tools](https://github.com/broadinstitute/warp-tools/tree/develop/3rd-party-tools/scvi-scanvi), whose changelog lists each image version and its key libraries (scvi-tools 1.5.1 and snapatac2 2.10 are pinned; the rest float within them). Every output records the versions that wrote it in `uns['scanvi_versions']`.

## Versioning

All scANVI pipeline releases are documented in the [scANVI changelog](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/scanvi/scANVI.changelog.md).

## Citing the scANVI Pipeline

If you use the scANVI Pipeline in your research, please identify the pipeline in your methods section using the [scANVI SciCrunch resource identifier](https://rrid.site/resolver/SCR_028705).

- Ex: *scANVI Pipeline (RRID:SCR_028705)*

When citing WARP, please use the following:

Kylee Degatano, Aseel Awdeh, Robert Sidney Cox III, Wes Dingman, George Grant, Farzaneh Khajouei, Elizabeth Kiernan, Kishori Konwar, Kaylee L Mathews, Kevin Palis, Nikelle Petrillo, Geraldine Van der Auwera, Chengchen (Rex) Wang, Jessica Way. "Warp Analysis Research Pipelines: Cloud-optimized workflows for biological data processing and reproducible analysis." _Bioinformatics_, 2025; [https://doi.org/10.1093/bioinformatics/btaf494](https://doi.org/10.1093/bioinformatics/btaf494)

Please also cite the underlying scvi-tools models, and scArches when samples were mapped onto a reference-only model:

- Lopez, R., Regier, J., Cole, M.B., Jordan, M.I., Yosef, N. "Deep generative modeling for single-cell transcriptomics." _Nature Methods_ 15, 1053–1058 (2018).
- Xu, C., Lopez, R., Mehlman, E., Regier, J., Jordan, M.I., Yosef, N. "Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models." _Molecular Systems Biology_ 17, e9620 (2021).
- Lotfollahi, M., Naghipourfar, M., Luecken, M.D., et al. "Mapping single-cell data to reference atlases by transfer learning." _Nature Biotechnology_ 40, 121–130 (2022).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

## Annotated h5ad spec

The format of each `<input file name>_annotated.h5ad`: one per GEX input, and one per ATAC input in multiome mode, compressed with gzip. n is the number of cells.

| Slot | Type | Contents |
| --- | --- | --- |
| `obs_names` | str | Input cell barcodes: GEX cells called by STARsolo (`star_IsCell`) with at least 3 counts; in multiome mode, only barcodes present in both GEX and ATAC. ATAC rows use the matching GEX barcode. |
| `var_names` | str | GEX: input genes with at least 3 counts. ATAC: the genes of the genome annotation. |
| `.X` | CSR, float32 | log2(CPM + 1) |
| `.raw.X` | CSR, float32 | Counts (ATAC: gene-activity counts) |
| `obs`, input columns | as input | Kept unchanged |
| `obs['n_counts']`, `var['n_counts']` | float32 | GEX only: counts per cell and per gene |
| `obs['batch']`, `obs['modality']` | category | Tags used in training: `pd-multiome_sci_gex` and `rna_unannotated` (ATAC: `pd-multiome_sci_atac` and `atac_unannotated`) |
| `obs['<level>']` | category | Predicted label at each level, for example `obs['subclass']` |
| `obs['<level>_max_probability']` | float32 | Maximum SCANVI posterior; omitted when `output_max_probability` is `false` |
| `obsm['X_scVI_<level>']` | float32, n × 30 | scVI latent at each level |
| `obsm['X_umap_<level>']` | float32, n × 2 | UMAP of that latent, fitted on this sample's cells |
| `obsm['X_umap']` | float32, n × 2 | Copy of the first level's UMAP |
| `uns`, input keys | as input | Kept, for example Optimus's `input_id` and `pipeline_version` |
| `uns['log1p']` | dict | `{'base': 2}` |
| `uns['scanvi_input_id']`, `uns['scanvi_input_file']`, `uns['scanvi_model']` | str | Sample id, GEX file name, model name |
| `uns['scanvi_levels']` | str array | Levels, in order |
| `uns['scanvi_versions']` | dict | Library versions of the image: scvi-tools, torch, scanpy, anndata, snapatac2, numpy, scipy |
| Not written | | `layers`, `obsp`, `varm`, `varp`, sample metadata |
