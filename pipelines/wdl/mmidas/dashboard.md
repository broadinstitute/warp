### MMIDAS: Mixture Model Inference with Discrete-coupled AutoencoderS

This workspace contains a set of example workflows that run **MMIDAS**, an unsupervised method for discovering reproducible cell types (and their continuous within-type variation) from single-cell/single-nucleus datasets. MMIDAS combines a generalized mixture model with a multi-armed deep neural network to *jointly* infer a discrete cell-type category and a continuous, type-specific variability ("state") for each cell. In this implementation, coupled mixture variational autoencoders (cpl-mixVAE) keep only the categories that independent encoder "arms" agree on, then iteratively prune categories until the surviving set is highly reproducible. The result is a set of discrete, consensus cell-type categories plus a continuous state that captures variation *within* each type.

MMIDAS was developed and published by **Yeganeh Marghi, Rohan Gala, Fahimeh Baftizadeh, and Uygar Sümbül** (Allen Institute for Brain Science). In their paper they show that modeling variability as continuous latent factors followed by a separate clustering step — or clustering the data directly — can make qualitative mistakes when the number of cell types is very large (hundreds to thousands), and they demonstrate MMIDAS on four brain single-cell datasets spanning different technologies, species, and conditions, in both unimodal and multimodal settings. The method and all scientific credit belong to those authors; please see [Citation and Credit](#citation-and-credit) and cite their work if you use these workflows.

The workflows are provided as a worked, end-to-end example on a public reference dataset. They are designed to run in the order presented, but each workflow can be launched independently using the sample inputs provided in this workspace.

The WDL workflow wrappers and this workspace were developed by the Data Sciences Platform at the Broad Institute; the underlying MMIDAS method and algorithms are the work of the original authors cited below.

> **Read this first — the data-prep workflow is an example, not a general-purpose tool.**
> `MMIDAS_DataPrep` is a *reference ingest* written specifically for the Allen Brain Atlas Smart-seq Mouse ALM/VISp files used here. **You will not be able to run your own dataset through it as-is.** The two model workflows (`MMIDAS_Train` and `MMIDAS_Analyze`) *are* general: they run on any AnnData `.h5ad` file that follows the simple contract described in [Bringing your own data](#bringing-your-own-data-replacing-mmidas_dataprep). Treat `MMIDAS_DataPrep` as a template you replace, and the other two workflows as the reusable engine.

---

## Workflows Overview

<!-- Optional: add an overview diagram here, e.g.
![](https://storage.googleapis.com/.../MMIDAS_overview.png)
-->

This workspace has three example workflows:

1. **MMIDAS_DataPrep** *(example ingest — dataset-specific, not general)*: converts the raw Allen Brain Atlas Smart-seq exon-count CSVs into a single normalized, filtered AnnData `.h5ad` file ready for training. This is the only dataset-specific stage; see [Bringing your own data](#bringing-your-own-data-replacing-mmidas_dataprep).

2. **MMIDAS_Train**: optionally trains a data augmenter, trains the core cpl-mixVAE model with iterative category pruning, and evaluates all checkpoints to recommend the optimal number of categories (`model_order`). It stops at a **human-review checkpoint**: you inspect the evaluation results and figures before continuing.

3. **MMIDAS_Analyze**: takes the reviewed model from `MMIDAS_Train` and produces the downstream biology — classification/clusterability analysis (how separable the discovered categories are) and state-traversal figures (what varies continuously within each category).

```
raw CSVs ──► MMIDAS_DataPrep ──► .h5ad ──► MMIDAS_Train ──► (human review) ──► MMIDAS_Analyze ──► figures
             (example only)                 model + eval                        classification +
                                            checkpoints                          state traversal
                                                  │                                    │
                                                  └────────────┬───────────────────────┘
                                                               ▼
                                              MMIDAS_output_validation.ipynb
                                              (checks the whole chain end to end)
```

These workflows cover MMIDAS's **transcriptomic** analysis, matching the authors' reference
notebooks. MMIDAS also supports coupled multimodal (expression + electrophysiology) analysis, which
is not implemented here — see
[Multimodal analysis](#multimodal-analysis-transcriptomics--electrophysiology).

A notebook, `MMIDAS_output_validation.ipynb`, validates a completed set of runs — see
[Validating a run](#validating-a-run--mmidas_output_validationipynb). Before your first run on your
own data, read
[Tuning for your own data](#tuning-for-your-own-data--read-this-before-your-first-run): three of the
training defaults are specific to the example dataset and will silently produce a useless model if
carried over unchanged.

---

## Sample Data

The example data is the **Allen Brain Atlas 2018 Mouse Smart-seq** dataset covering two cortical regions — primary visual cortex (**VISp**) and anterior lateral motor cortex (**ALM**). It consists of full-length Smart-seq exon-count matrices plus per-cell metadata (including reference cell-type "cluster" labels curated by the Allen Institute).

The raw files expected by `MMIDAS_DataPrep` are:

| File | Description |
| --- | --- |
| `mouse_VISp_2018-06-14_exon-matrix.csv` | VISp raw exon counts (genes × cells) |
| `mouse_VISp_2018-06-14_samples-columns.csv` | VISp per-cell metadata (one row per cell, same order as matrix columns) |
| `mouse_ALM_2018-06-14_exon-matrix.csv` | ALM raw exon counts (genes × cells) |
| `mouse_ALM_2018-06-14_samples-columns.csv` | ALM per-cell metadata |
| `mouse_ALM_2018-06-14_genes-rows.csv` | Full gene list (must contain a `gene_symbol` column) |
| `genes_SS_ALM-VISp.csv` | Selected 5,032-gene subset used for training (must contain a `genes` column) |

Optional reference files used by `MMIDAS_Analyze`:

| File | Used for |
| --- | --- |
| `tree_Mouse_ALM-VISp_2018.csv` | Hierarchical taxonomy tree — orders categories by their dominant reference cell type (optional) |
| `KEGG.toml` | KEGG pathway gene sets — enables per-pathway box plots in the state-traversal step (optional) |

These files are provided in the workspace bucket and referenced by the example input JSONs.

---

## Workflows

### 1. MMIDAS_DataPrep  *(example ingest — dataset-specific)*

**What does it do?**

`MMIDAS_DataPrep` runs `01_data_prep.py`, which turns the raw Allen Smart-seq CSVs into one analysis-ready AnnData `.h5ad` file. In detail it:

1. Loads the VISp and ALM exon-count matrices, reading **only** neuronal-cell columns to keep memory low.
2. Retains only the requested neuronal classes (default `GABAergic` and `Glutamatergic`).
3. Concatenates the two regions into a single matrix.
4. Normalizes counts to **log-CPM**: `log1p(counts / rowsum × 1e6)`.
5. Subsets to the selected gene list.
6. Removes low-quality / rare clusters (default `Low Quality,CR Lhx5,Meis2 Adamts19`).
7. Applies two dataset-specific cell-type renames to match the reference taxonomy.
8. Writes the result as a `.h5ad` with the expression matrix in `X`, cell metadata in `obs` (including the reference `cluster` labels), and gene symbols as `var_names`.

> **Why this is an example and not a general tool.** Almost every step above encodes assumptions that are specific to this dataset: exactly two regions concatenated together; a fixed CSV layout with *positional* alignment between the count matrix columns and the metadata rows; hard-coded column names (`class`, `cluster`, `gene_symbol`, `genes`); a fixed log-CPM normalization that assumes raw-count input; and hard-coded cluster-removal and rename lists. **Your own data will almost certainly have a different raw format, different metadata columns, and different QC choices, so it cannot flow through this workflow unchanged.** This is expected — data ingest is inherently dataset-specific. See [Bringing your own data](#bringing-your-own-data-replacing-mmidas_dataprep) for the output contract you need to reproduce.

**What does it require as input?**

| Input | Type | Description |
| --- | --- | --- |
| `visp_exon_matrix` | File | VISp raw exon count matrix CSV (genes × cells) |
| `visp_samples` | File | VISp per-cell metadata CSV |
| `alm_exon_matrix` | File | ALM raw exon count matrix CSV |
| `alm_samples` | File | ALM per-cell metadata CSV |
| `genes_rows` | File | Full gene list CSV (with `gene_symbol` column) |
| `selected_genes` | File | Selected gene subset CSV (with `genes` column) |
| `output_basename` | String | Basename for the output `.h5ad` (default `Mouse_ALM-VISp_cpm`) |
| `remove_clusters` | String | Comma-separated clusters to exclude (default `Low Quality,CR Lhx5,Meis2 Adamts19`) |
| `neuronal_classes` | String | Comma-separated cell classes to keep (default `GABAergic,Glutamatergic`) |

**What does it return as output?**

| Output | Type | Description |
| --- | --- | --- |
| `preprocessed_h5ad` | File | The normalized, filtered AnnData `.h5ad` — the input to `MMIDAS_Train` |
| `pipeline_version_out` | String | Pipeline version string |

---

### 2. MMIDAS_Train

**What does it do?**

> **This workspace is configured for a shorter run than the published study.**
> `n_epoch_p` (epochs per pruning round) is set to **1,000** here, giving 52,000 total epochs —
> **~13.5 hours, ~$13**. The published analysis used **10,000**, which is 430,000 epochs —
> **~4.5 days, ~$110**. The example outputs in this workspace, and the validation-notebook results
> reported below, all come from the 1,000-epoch configuration; on the example data it reached
> `avg_consensus 0.969` and `model_order 89` against the published 92.
>
> **To run the published configuration**, set `n_epoch_p = 10000` in the workflow inputs (or use
> `example_inputs/MMIDAS_Train.json`, which keeps that value) and plan for the longer runtime and
> cost. `max_prun_it` is 42 in both, so the published `model_order` is reachable either way.

`MMIDAS_Train` is the core modeling stage. It runs three steps and ends at a human-review checkpoint:

- **(Optional) Augmenter training** — trains a UDAGAN VAE-GAN augmenter that can generate realistic synthetic cells to stabilize training. Off by default (`run_augmenter = false`).
- **cpl-mixVAE training with pruning** — trains two coupled encoder arms starting from an upper bound of `n_categories` categories, then iteratively prunes the least-reproducible category (lowest inter-arm consensus) for up to `max_prun_it` rounds. Each round writes a model checkpoint.
- **Evaluation** — scores every checkpoint, runs K-selection to recommend the optimal number of categories (`model_order`), and writes `evaluation_results.json` plus consensus/K-selection figures.

**Human-review checkpoint:** after this workflow finishes, download `evaluation_results.json` and the evaluation figures and confirm the recommended `model_order` is biologically sensible **before** launching `MMIDAS_Analyze`. The evaluation JSON and the model tarball are the hand-off files to the next stage.

Three fields in `evaluation_results.json` decide whether the run is usable at all, and none of them is `model_order`:

| Field | Reject the run if |
| --- | --- |
| `k_selection_met_threshold` | `false` — no checkpoint reached `k_select_thr`, so `model_order` came from a fallback rather than a selection |
| `n_populated_categories` | far below `model_order` — `model_order` counts categories that survived pruning, which stays high even when the model routes every cell into a handful of them |
| `collapse_warning` | non-null — the two above disagree badly enough that downstream Analyze figures will be dominated by empty categories |

`avg_consensus` should also be at or above `k_select_thr`. A run with high `model_order` and near-zero `avg_consensus` has not found reproducible categories; it has failed to train its discrete latent. In the training log, watch the per-epoch `Entropy` against the `uniform=` value printed next to it — an `Entropy` that stays pinned at `uniform` while the reconstruction loss falls means the categorical variable never committed and no amount of pruning will fix it. The usual cause is `tau` being too large for your `n_categories`; see [Tuning for your own data](#tuning-for-your-own-data--read-this-before-your-first-run).

After `MMIDAS_Analyze` completes, run `MMIDAS_output_validation.ipynb` to check the whole chain automatically rather than eyeballing the JSON — see [Validating a run](#validating-a-run--mmidas_output_validationipynb).

**Detail for the detail-inclined.** The model is a *coupled* mixture VAE: two (or more) arms encode the same cell independently, and the training loss penalizes disagreement between the arms' categorical assignments (`lam`/`lam_pc` coupling factors). Only categories the arms agree on survive pruning, which is what makes the discovered categories reproducible rather than an artifact of a single run — this consensus-across-arms idea is the core contribution of the MMIDAS method (Marghi et al., 2024; see [Citation and Credit](#citation-and-credit)). Each cell also gets a low-dimensional continuous **state** variable (`state_dim`) that captures within-type variation. Reconstruction can use `MSE` or `ZINB` loss (`training_mode`).

**What does it require as input?**

Key inputs (all hyperparameters have production defaults):

| Input | Type | Default | Published value | Description |
| --- | --- | --- | --- | --- |
| `preprocessed_h5ad` | File | — | — | Output of `MMIDAS_DataPrep`, or your own contract-compliant `.h5ad` |
| `run_augmenter` | Boolean | `false` | `false` | Whether to train the optional data augmenter first |
| `n_categories` | Int | `120` | `120` | Upper-bound number of categories before pruning |
| `n_arm` | Int | `2` | `2` | Number of coupled encoder arms |
| `state_dim` | Int | `2` | `2` | Continuous within-type state dimension |
| `latent_dim` | Int | `10` | `10` | Low-dimensional embedding dimension |
| `training_mode` | String | `MSE` | `MSE` | Reconstruction loss (`MSE` or `ZINB`) |
| `n_epoch` | Int | `10000` | `10000` | Epochs before pruning begins |
| `n_epoch_p` | Int | **`1000`** | **`10000`** | Epochs per pruning round. **The one default that differs from the published study** — see the note above |
| `max_prun_it` | Int | `42` | `42` | Maximum pruning iterations |
| `min_con` | Float | `0.99` | `0.99` | Reporting only in this implementation; pruning runs the full `max_prun_it` regardless |
| `tau` | Float | `0.005` | `0.005` | Categorical softmax temperature. **Rescale this if you change `n_categories`** — see [Tuning for your own data](#tuning-for-your-own-data--read-this-before-your-first-run) |
| `k_select_thr` | Float | `0.95` | `0.95` | Consensus threshold used to recommend `model_order` |
| `batch_size` | Int | `5000` | `5000` | Mini-batch size |
| `seed` | Int | `0` | unseeded | Random seed; the reference train/test split was unseeded |
| `train_gpu` | Int | `0` | — | Set to `1` to attach a GPU (see note below) |

> **GPU note.** Training is much faster on a GPU, and every figure in [Time and Cost Estimates](#time-and-cost-estimates) assumes one. Set `train_gpu = 1` and that is all: the `TrainMixVAE` task's runtime block already declares `gpuCount` and `gpuType: "nvidia-tesla-t4"` and passes `--cuda` to the training script. There is **no GPU setting to enable in Terra** for workflow submissions — Terra's Cloud Environment has GPU options, but those apply to interactive notebooks and RStudio, not to Cromwell workflow tasks.
>
> The one thing that can block a GPU task is **GCP quota**: the Google project behind your Terra billing project needs available GPU quota in the execution region, or the task will fail to schedule rather than fall back to CPU. If that happens, raise it with Terra support — billing-project quotas are not adjustable from the Terra UI.

**What does it return as output?**

| Output | Type | Description |
| --- | --- | --- |
| `evaluation_results_json` | File | **Review this.** Recommended `model_order`, selected model, and metrics |
| `evaluation_figures` | Array[File] | Consensus heatmaps and K-selection curves for review |
| `summary_performance` | File | Per-checkpoint consensus / reconstruction pickle behind the K-selection decision |
| `checkpoints_manifest` | File | Manifest of all model checkpoints and architecture settings |
| `model_tar` | File | Tarball of all trained checkpoints |
| `augmenter_checkpoint` | File? | Augmenter model (only if `run_augmenter = true`) |
| `pipeline_version_out` | String | Pipeline version string |

---

### 3. MMIDAS_Analyze

**What does it do?**

`MMIDAS_Analyze` takes the reviewed model from `MMIDAS_Train` and produces the downstream biological analyses. It first restores the model checkpoints, then runs two analyses in parallel and turns each into figures:

- **Clusterability (steps 03b → 04)** — trains a random-forest classifier and computes silhouette scores to quantify how separable the MMIDAS categories are, compared against a PCA baseline and the reference cluster labels. Produces classification-accuracy bar charts, silhouette curves, and confusion-matrix heatmaps.
- **State traversal (steps 03c → 05)** — walks along each category's continuous state axis and visualizes how gene expression changes, producing per-category state-space scatter plots and (if a KEGG file is supplied) per-pathway box plots.

**Detail for the detail-inclined.** The "clusterability" analysis answers *"are these categories real and separable?"* by asking how well a classifier can recover them and how tight/separated they are in embedding space. The "state traversal" answers *"what varies continuously within a type?"* by holding the categorical assignment fixed and moving along the state latent, then decoding the resulting expression profiles. Category ordering can optionally follow the Allen hierarchical taxonomy tree (`htree_file`).

**What does it require as input?**

| Input | Type | Description |
| --- | --- | --- |
| `preprocessed_h5ad` | File | Same `.h5ad` used for training |
| `checkpoints_manifest` | File | From `MMIDAS_Train` |
| `model_tar` | File | From `MMIDAS_Train` |
| `evaluation_results_json` | File | From `MMIDAS_Train` — **after** you have reviewed it |
| `kegg_toml` | File? | Optional KEGG pathway file (enables pathway box plots) |
| `htree_file` | File? | Optional taxonomy tree (enables taxonomy-ordered categories) |
| `n_pca` | Int | PCA components for the linear baseline (default `100`) |
| `k_fold` | Int | Cross-validation folds for classification (default `10`) |
| `n_traversal_steps` | Int | Points along each state traversal (default `50`) |
| `traversal_arm` | Int | Which encoder arm to visualize (default `0`) |
| `n_selected_cats` | Int | Number of categories to plot, `0` = all (default `10`) |
| `batch_size` / `seed` | Int | Must match training (defaults `5000` / `0`) |

**What does it return as output?**

| Output | Type | Description |
| --- | --- | --- |
| `clusterability_figures` | Array[File] | Classification accuracy, silhouette, and confusion-matrix figures |
| `clusterability_manifest` | File | Manifest for the clusterability outputs |
| `state_traversal_figures` | Array[File] | Per-category state-space and (optional) pathway figures |
| `state_traversal_manifest` | File | Manifest for the state-traversal outputs |
| `pipeline_version_out` | String | Pipeline version string |

---

## Bringing Your Own Data (replacing MMIDAS_DataPrep)

You cannot run your own dataset through `MMIDAS_DataPrep` — it is hard-wired to the Allen ALM/VISp CSV format. Instead, produce your own AnnData `.h5ad` with any tool you like (e.g. Scanpy) that satisfies the small contract below, then feed it straight into `MMIDAS_Train` (and pass the same file to `MMIDAS_Analyze`).

**The `.h5ad` contract expected by MMIDAS_Train / MMIDAS_Analyze:**

| Where | Requirement |
| --- | --- |
| `adata.X` | Log-normalized expression matrix, cells × genes (the example uses log-CPM: `log1p(counts / rowsum × 1e6)`). A sparse `float32` matrix is recommended. |
| `adata.var_names` | Gene symbols/identifiers, one per column of `X`. |
| `adata.obs['cluster']` | A **reference cell-type label per cell** (string). This is used as the ground-truth label for evaluation and classification. This column must exist. |
| `adata.obs['subclass']`, `adata.obs['class']` | Optional additional label columns used for some ordering/plots. |

Notes:

- The model itself is dataset-agnostic. Everything that would differ per dataset — number of genes, number of categories (`n_categories`), embedding sizes, epochs — is a workflow parameter, so you tune those to your data rather than editing code.
- The reference `cluster` labels are used to *evaluate* and *order* the discovered categories; they are not required for the model to learn, but the evaluation, classification, and taxonomy-ordering steps expect them.
- `01_data_prep.py` in the [warp-tools](https://github.com/broadinstitute/warp-tools) repo is a good worked example of how to build a contract-compliant `.h5ad`; copy and adapt it for your own raw format.
- `MMIDAS_output_validation.ipynb` validates against the authors' published Mouse ALM/VISp results, so its reference comparisons will not apply to your data. Adapt it as a starting point for your own validation — see [Using it on your own data](#using-it-on-your-own-data).

---

## Multimodal analysis (transcriptomics + electrophysiology)

MMIDAS as published is not limited to gene expression. The authors describe it as applying to "both,
uni-modal and multi-modal datasets," and the paper demonstrates coupled analysis of **transcriptomic
and electrophysiological** measurements from the same cells — Patch-seq data, where each neuron is
both patched for its electrical properties and sequenced. The appeal is that consensus is then
required *across modalities*: a cell type is kept only if it is recoverable from both what a neuron
does electrically and what it expresses, which is a stronger claim than either modality alone
supports.

> **These three workflows implement the transcriptomic path only.** They correspond to the authors'
> reference notebooks (`1_data_prep` … `5_state_traversal`), which are single-modality. Nothing in
> this workspace has been run against electrophysiology data, and none of it is tested for that. The
> section below is a pointer for anyone who wants to go that direction, not a supported path.

**Why it is not just a matter of different inputs.** The multimodal model is a different model, not
the same one with an extra file:

- The `.h5ad` contract above describes a single expression matrix. A multimodal run needs two feature
  matrices plus the pairing that says which row of each belongs to the same cell.
- `MMIDAS_Train` calls `cpl_mixVAE.init_model()`, whose signature takes `n_arm` — the number of
  encoder arms over one modality — and has no modality parameter. The network class it builds is
  `RNA_RNA_mixVAE`, i.e. arms coupled across one data type.
- The multimodal model instead needs per-modality architecture and coupling: separate arm counts and
  state dimensions for each modality, and a cross-modal coupling factor in addition to the
  within-modality ones.

**Where to start.** The authors' repository includes `tutorials/datasets_training/train_patchseq.py`,
which exposes exactly those parameters — `--n_modal`, `--n_arm_T` / `--n_arm_E`, `--state_dim_T` /
`--state_dim_E`, and `--lam_T` / `--lam_E` / `--lam_TE` for transcriptomic, electrophysiological, and
cross-modal coupling. Read it as the specification for what a multimodal workflow would need to
provide. Note that it imports helper modules (`utils.training`, `utils.helpers`, including a
`load_patchseq` loader) that are not part of the revision pinned in this workspace's Docker image, so
those pieces would need to be sourced before it will run.

**Data and feature extraction.** The authors use the
[Allen Institute Patch-seq dataset](https://dandiarchive.org/dandiset/000020/) and note that the
electrophysiological features were computed following the approach in their companion `cplAE_MET`
repository. Deriving those features from raw traces is a substantial preprocessing step in its own
right, and is the analogue of `MMIDAS_DataPrep` for the electrophysiology side.

**Scope of the work required.** Realistically this is a fourth workflow rather than an option on the
existing one: an ingest stage producing paired transcriptomic and electrophysiological matrices, a
training stage wrapping the multimodal model, and evaluation and analysis stages that report
consensus per modality as well as across them. `MMIDAS_Analyze` would need corresponding changes,
since its classification and state-traversal steps assume a single feature space.

---

## Tuning for your own data — read this before your first run

The defaults in these workflows are tuned for the example dataset (Mouse ALM/VISp, 22,365 cells,
5,032 genes, 115 reference t-types, `n_categories = 120`). Three of them are **not** safe to carry
over to a different dataset or a different `n_categories`, and getting them wrong produces a run
that completes successfully and reports plausible-looking numbers while being useless.

### 1. `tau` must be rescaled whenever you change `n_categories`

This is the single most important parameter to get right. `tau` is the categorical softmax
temperature, and `cpl_mixVAE.init_model` documents it as *"usually equals to 1/n_categories"*.

| `n_categories` | Appropriate `tau` |
| --- | --- |
| 15 | ~0.067 |
| 50 | ~0.020 |
| **120 (this workspace)** | **~0.008 — default is 0.005** |
| 250 | ~0.004 |

If `tau` is too large for your `n_categories`, the categorical posterior stays nearly flat and the
model reconstructs entirely through the continuous state variable rather than the discrete categories.
The effect is large: on the example data at `n_categories = 120`, a `tau` sized for
`n_categories = 15` gave `avg_consensus` **0.026** and **11** populated categories, where the correctly
scaled `tau = 0.005` gave `avg_consensus` **0.969** and 89 of 89 populated. The failure is silent —
nothing errors, and `model_order` still comes back a plausible number.

**How to tell within the first few hours.** The `TrainMixVAE` log prints, every epoch:

```
Entropy: -0.8030 (uniform=-9.5750)
```

`uniform` is the value `Entropy` takes when both arms' categorical posteriors are completely flat —
i.e. the discrete latent carries no information. Watch the gap:

- `Entropy` **moving decisively away from `uniform`** → the categorical variable is committing. Good.
- `Entropy` **pinned near `uniform`** while the reconstruction loss falls → collapse. Kill the run
  and lower `tau`. This is visible at the end of the pre-pruning phase, roughly 2.5 hours in, long
  before pruning starts.

For reference, the healthy run moved from −9.57 to −0.80 (about 1.5 effective categories per cell);
the collapsed run only reached −9.02 (about 91 of 120 — essentially no commitment).

This is the cheapest possible check on a long run: it is readable at the end of the pre-pruning phase
(the first `n_epoch = 10000` epochs, ~2.5 h, ~$2.50), before any pruning starts. If `Entropy` is
still pinned near `uniform` there, lower `tau` and restart rather than paying for 42 pruning rounds.

### 2. `max_prun_it` bounds which answers are even reachable

Pruning removes **one** category per round, so the smallest `model_order` a run can produce is:

```
n_categories - max_prun_it
```

If the number of cell types in your data falls below that floor, no amount of training will find it —
the answer is outside the search space. With the defaults (`n_categories = 120`,
`max_prun_it = 42`) the reachable range is **78–120**. If you expect ~30 types from
`n_categories = 120`, you need `max_prun_it` ≥ 90.

Set `n_categories` generously above your expected type count and `max_prun_it` large enough that
your plausible range sits comfortably inside the reachable window.

### 3. Runtime and cost scale with `max_prun_it × n_epoch_p`

Total epochs are `n_epoch + max_prun_it × n_epoch_p`, and training dominates everything else in the
pipeline. Both parameters therefore multiply your bill: raising `max_prun_it` to widen the search
space (above) also raises the cost proportionally. See
[Time and Cost Estimates](#mmidas_train) for the measured figures.

Two things worth knowing before you raise `n_epoch_p` from its default of 1,000 to the published
10,000: on the example data the 1,000-epoch configuration already reached `avg_consensus 0.969` and
`model_order 89` against the published 92, and round-by-round consensus in a full-length run
plateaued by round 2 and then moved only within noise for 17 more rounds.

### 4. There is no resume — protect against losing a long run

`TrainMixVAE` writes `model.tar.gz` only when the task **completes**. If it is aborted (a cost cap,
a timeout, a manual cancel), the intermediate checkpoints are lost with the VM and the run must
start over. `preemptible` is already `0` so the VM will not be reclaimed mid-run, but:

- **Set Terra cost caps above the expected spend** — ≥ ~$25 at the default `n_epoch_p = 1000`,
  ≥ ~$150 at the published 10,000. A cap hit mid-run discards the whole run, not just the remainder.
- 4–5 days is within the usual 7-day GCP task ceiling, but only just. Confirm your project does not
  impose a shorter limit before launching at `n_epoch_p = 10000`.

### 5. Smaller things worth knowing

| Parameter / behaviour | What to know |
| --- | --- |
| `min_con` | **Reporting only** in this implementation: pruning runs the full `max_prun_it` regardless of the value set. Do not expect `min_con` to halt anything. |
| `k_select_thr` | If no checkpoint reaches it, `K_selection` returns nothing and `Evaluate` falls back to the un-pruned checkpoint. Check `k_selection_met_threshold` in `evaluation_results.json` — `false` means `model_order` came from a fallback, not a selection. |
| `kegg_toml` | Optional. Omit it and `n_pathways` is 0 with no pathway figures — expected, not a failure. Supply it and confirm `n_pathways > 0`; zero pathways *with* a `kegg_toml` means gene-name lookup failed. |
| `htree_file` | Optional; enables taxonomy ordering in stage 03c. |
| `n_selected_cats` | Capped at the number of *populated* categories, so the manifest may report fewer than you asked for. |
| Run-to-run variation | Training is **not** bit-reproducible even with a fixed `seed`. Two runs of the identical configuration gave `model_order` 96 and 89 with `avg_consensus` 0.9075 and 0.9692. Do not build a conclusion on one run — see [Run-to-run variation](#run-to-run-variation). |
| `Classify` retries | This task has retried on three consecutive runs (10-fold random forest over the full cell set). It succeeds on retry, but its outputs land in `call-Classify/attempt-N/` rather than `call-Classify/`. |

---

## Running the Workflows

The workflows are pre-configured with the example inputs in this workspace (see the `example_inputs/` JSON files). For each workflow:

1. Select the workflow from the **Workflows** tab.
2. Provide inputs — either use the provided example JSON or edit the input fields.
3. (For `MMIDAS_Train`) enable a GPU by setting `train_gpu = 1` for a much faster run.
4. Launch the workflow.

Recommended order:

1. Run **MMIDAS_Train** on the provided example `.h5ad` (or your own contract-compliant `.h5ad`).
2. **Review** `evaluation_results.json` and the evaluation figures; check `k_selection_met_threshold`, `n_populated_categories` and `collapse_warning`, then confirm `model_order`.
3. Run **MMIDAS_Analyze**, passing the `checkpoints_manifest`, `model_tar`, and reviewed `evaluation_results_json` from step 1.

If you want to reproduce the example end-to-end from the raw Allen CSVs, run **MMIDAS_DataPrep** first to produce the `.h5ad` — but remember this step only works for the Allen ALM/VISp files.

---

## Time and Cost Estimates

Measured on Terra with the example dataset (22,365 cells × 5,032 genes). These are the figures every
other section of this document refers back to. Training dominates; the other two stages are minor by
comparison.

| Stage | Configuration | Time | Cost |
| --- | --- | --- | --- |
| `MMIDAS_DataPrep` | — | ~9 min | < $1 |
| **`MMIDAS_Train`** | **`n_epoch_p = 1000`** (workspace default, 52,000 epochs) | **~13.5 h** | **~$13** |
| `MMIDAS_Train` | `n_epoch_p = 10000` (published, 430,000 epochs) | ~4.5 days | ~$110 |
| `MMIDAS_Analyze` | — | ~1.4 h | < $2 |
| | **End-to-end at the default** | **~15 h** | **~$15** |

### MMIDAS_Train

Training cost is set almost entirely by total epochs, `n_epoch + max_prun_it × n_epoch_p`, at
**~0.9 s/epoch** on an `nvidia-tesla-t4` at **~$1.00/hr**. Those two constants are what the table
above is derived from, so you can price any configuration from them.

The workspace default (`n_epoch_p = 1000`) and the published configuration (`10000`) differ in that
one field; `max_prun_it` is 42 in both. CPU-only training was not benchmarked and is impractical at
these epoch counts.

### MMIDAS_DataPrep

Needs `mem_size = 48` GiB: it holds the full 22,439 × 45,768 count matrix in memory before subsetting
to the selected genes.

### MMIDAS_Analyze

CPU only, no GPU; the `Classify` task accounts for most of the 1.4 hours. Cheap enough to re-run
freely, which matters because it is the stage you re-run when figures or downstream analysis change
without retraining.

For more information about controlling Cloud costs, see [this article](https://support.terra.bio/hc/en-us/articles/360029748111).

---

## Fidelity to the original MMIDAS analysis

These workflows are a port of the original MMIDAS code, recreating the authors' notebook analysis as
WDL workflows. The [MMIDAS repo](https://github.com/AllenInstitute/MMIDAS) ships its notebooks with
outputs saved, so the authors' results for the example dataset are recorded and can be compared
against directly.

**Reference values for Mouse ALM/VISp:**

| Quantity | Reference | Source |
| --- | --- | --- |
| matrix shape | 22,365 × 5,032 | `1_data_prep.ipynb` |
| reference t-types | 115 | `2_train.ipynb` data summary |
| pruning rounds | 42 | `3_evaluation.ipynb` (checkpoints `after_pruning_1..42`) |
| `model_order` | **92** | `3_evaluation.ipynb`; hardcoded in notebooks 4 and 5 |
| `avg_consensus` | 0.939 (test cells) / 0.954 (K-selection) | `3_evaluation.ipynb` |

**Measured results in this workspace.** Two independent runs of the *identical* configuration are
shown, because the difference between them is itself a result you need to know about — see
[Run-to-run variation](#run-to-run-variation) below.

| Quantity | Reference | Run A | Run B | |
| --- | --- | --- | --- | --- |
| matrix shape | 22,365 × 5,032 | 22,365 × 5,032 | 22,365 × 5,032 | match |
| reference t-types | 115 | 115 | 115 | match |
| pruning rounds | 42 | 42 | 42 | match |
| `model_order` | 92 | **96** | **89** | within tolerance either way |
| `avg_consensus` | 0.939 / 0.954 | **0.9075** | **0.9692** | reference falls between them |
| populated categories | 92 (all) | **96 of 96** | **89 of 89** | all populated, both arms |
| K-selection met `k_select_thr` | yes | **yes** | **yes** | no fallback |

Run A is the run distributed with this workspace. Run B was produced earlier in a development
workspace and is included here only as the second sample.

`MMIDAS_output_validation.ipynb` checks a run against these under **Stage 6 — Reference
comparison**; see [Validating a run](#validating-a-run--mmidas_output_validationipynb) for how to
point it at your own submissions and what the check labels mean. Where the results *do* differ, see
[Where the results differ from the published analysis](#where-the-results-differ-from-the-published-analysis).

### Run-to-run variation

**Two runs of the same configuration on the same input will not give the same answer.** Plan for
this before you build any conclusion on a single run.

| | `model_order` | `avg_consensus` |
| --- | --- | --- |
| Run A | 96 | 0.9075 |
| Run B | 89 | 0.9692 |
| Published | 92 | 0.939 / 0.954 |
| **Observed spread** | **7 categories** | **0.062** |

Both runs completed all 42 pruning rounds with every surviving category populated in both arms and
`k_selection_met_threshold` true, so both are valid runs — the spread is the algorithm's, not a
defect. The published values fall between the two, which is the reassuring part: the runs bracket the
reference rather than sitting to one side of it.

Two causes, neither removable:

- **The train/test split.** The reference calls `get_loaders` without a seed. These workflows seed it
  (`seed = 0`) so a given workflow run is at least self-consistent, but that does not recover the
  reference's split.
- **GPU non-determinism.** Floating-point reductions on a GPU are not associative, so identical
  inputs and an identical seed still diverge. In one pair of runs the two encoder arms swapped which
  one converged better.

Practical consequences:

- Do not treat a single `model_order` as *the* number of cell types. Report it as an estimate, and
  run the configuration more than once if the exact count matters to your conclusion.
- Thresholds in the validation notebook are set to admit this spread. `model_order_tol` is 5 pruning
  rounds and `avg_consensus_min` is 0.900, the latter calibrated from these two runs plus the
  published value. A collapsed model measures around 0.03, so the floor is nowhere near loose enough
  to pass a failed run.
- If two of your own runs differ by *much* more than the above — tens of categories, or consensus
  moving by several tenths — that is no longer ordinary variation. Check `tau` against your
  `n_categories` first.

### Where the training defaults come from

The reference notebooks pass only `n_categories`, `state_dim`, `n_arm` and `latent_dim` to
`cpl_mixVAE.init_model()` and let everything else fall to that function's defaults, so `MMIDAS_Train`
follows those defaults. The training-loop values (`n_epoch`, `n_epoch_p`, `max_prun_it`) come from
`tutorials/train_mixvae.py`, since `2_train.ipynb`'s `n_epoch = 10` is a walkthrough demo rather than
the configuration behind the published model.

Two are worth calling out:

- **`max_prun_it` bounds which answers are reachable.** Pruning removes one category per round, so
  the smallest `model_order` a run can produce is `n_categories - max_prun_it`. Reaching the
  reference `model_order` of 92 from `n_categories = 120` requires 28 rounds; the default is 42, the
  number the reference used.
- **`tau` scales with `n_categories`.** `init_model` documents it as "usually equals to
  1/n_categories" — about 0.0083 at `n_categories = 120`, and the default here is the library's
  0.005. If you change `n_categories`, rescale `tau` to match; see
  [Tuning for your own data](#tuning-for-your-own-data--read-this-before-your-first-run).

`n_epoch_p` is the one default that deliberately differs from the published configuration — 1,000
here against 10,000 — because it is what the workspace's example outputs were produced with. See the
note in [MMIDAS_Train](#2-mmidas_train).

`min_con` is **reporting only** in this implementation: pruning runs the full `max_prun_it` regardless
of its value, which is how every reference invocation behaves.

### Known divergences from the notebooks

Two differences are unavoidable in a generic workflow and are deliberate:

| Divergence | Why |
| --- | --- |
| **State-traversal category selection.** `5_state_traversal.ipynb` hardcodes `selected_c = [80, 119, 1, 13, 92, 25, 55, 110, 62, 69, 31]`. The workflow instead selects the `n_selected_cats` most-populated categories. | A workflow cannot reproduce a hand-picked list chosen by inspection. Selecting by population at least guarantees the figures show categories the model actually uses. |
| **Train/test split seeding.** `2_train.ipynb` calls `get_loaders` without a `seed`, so its split is random and unrecoverable. The workflow seeds it (`seed = 0`). | An unseeded split makes a workflow non-reproducible run to run. This is why an exact numerical match to the reference is not expected, and why the validation notebook compares `model_order` within a tolerance and `avg_consensus` against a range. |

### Where the results differ from the published analysis

The validated run reproduces the reference within tolerance. Three things do not match exactly, and
all three are worth knowing before you present results.

**1. `model_order` 96 against the published 92** — four pruning rounds apart, and 89 on the second
run. This is ordinary stochastic variation rather than a divergence in the port; see
[Run-to-run variation](#run-to-run-variation) for the spread and its causes. Exact agreement is not
achievable.

**2. t-type classification accuracy sits further below the PCA baseline than in the reference.**

| | PCA-100 | MMIDAS-10 | gap |
| --- | --- | --- | --- |
| Reference (`4_clusterability.ipynb`) | ~84.5% | ~73.5% | ~0.11 |
| Run A | 90.9% | 62.1% | 0.288 |
| Run B | 90.9% | 66.2% | 0.247 |

Every comparison points the same direction as the reference — PCA ahead on the 115 reference t-types,
the 10-dimensional MMIDAS embedding well ahead on MMIDAS's own categories (94.7% vs 64.8%) — so this
is a difference in magnitude on a downstream metric, not in behaviour. Some gap here is expected by
design: a 10-dimensional embedding is not attempting to beat a 100-component PCA basis at recovering
115 reference labels. Three caveats on the comparison itself: the metric has ~0.05 fold-to-fold
spread, it tracks `avg_consensus` across runs (the weaker-consensus run shows the wider gap), and
`4_clusterability.ipynb` reports these as a figure rather than numbers, so the reference column is
approximate. The validation notebook therefore reports this as *advisory*.

**3. The accuracy bar chart has fewer groups than the reference's.** The reference plots three label
sets — t-types, *Merged t-types*, and MMIDAS T Categories. The workflow plots t-types and the per-arm
MMIDAS categories, omitting the merged-t-type group. The two groups that do appear match the
reference's layout and direction.

---

## Validating a run — `MMIDAS_output_validation.ipynb`

The workspace includes a notebook that checks a completed set of runs end to end. Point it at your
three submissions, run it top to bottom, and it reports what passed, what failed, and which figures
still need a human eye.

**It validates against the authors' published results for the example dataset.** Its reference values
— matrix shape, 115 t-types, `model_order 92`, consensus, 42 pruning rounds, the accuracy figures —
are all read from the reference notebooks for Mouse ALM/VISp. That makes it a direct answer to "did
this run reproduce the published analysis?" on the example data, and it is the fastest way to know a
run is trustworthy before building analysis on top of it. On any other dataset those comparisons do
not apply; see [Using it on your own data](#using-it-on-your-own-data) below.

**Pointing it at your run.** The Config cell holds three constants — `DATAPREP_RUN`, `TRAIN_RUN`,
`ANALYZE_RUN` — each a Terra submission/workflow path. Everything else is derived, so repointing at
a new run is three lines. Figure lists are given as `gs://` prefixes and globbed, not enumerated.

**Two traps when copying paths from Terra:**

- A Terra output path contains **two** UUIDs — the submission ID and the workflow ID — and the bucket
  name contains a third (`fc-<uuid>`). Pasting the bucket's UUID where the workflow ID belongs
  produces a path that looks right and fails with an opaque `gsutil` error. If Stage 0 reports tasks
  as unreadable, check this first.
- `Classify` output lives under `call-Classify/attempt-N/` when that task retries, which it has done
  on every run so far. The notebook discovers the attempt automatically; you only supply the call
  directory.

**Checks are labelled by what they answer**, and only the first two decide the verdict:

| Kind | Question | Dataset-specific? | On failure |
| --- | --- | --- | --- |
| `plumbing` | Did the workflow execute correctly? Files present, shapes consistent, manifests mutually agreeing, all stages consuming the same inputs. | No — applies to any dataset | The port is broken. Fix before interpreting anything. |
| `fidelity` | Does this run reproduce the published analysis? Compared against values recorded in the reference notebooks. | **Yes — Mouse ALM/VISp only** | The run does not match how the authors ran it. |
| `advisory` | Is the model any good? Soft metrics with wide run-to-run spread, or comparisons against figures rather than published numbers. | Partly | Informational. Never fatal. |

On the example dataset, a run that reproduces the reference is a success even if advisory items
complain — and a run with better-looking numbers that does *not* reproduce the reference is not.

**Result for the run distributed with this workspace:** 47/48 — `plumbing 36/36`, `fidelity 7/7`,
`advisory 4/5`, verdict *"the workflow executed correctly and matches the reference analysis within
tolerance. The port is faithful."* The single advisory failure is the t-type accuracy gap described
above.

**What it cannot tell you.** The checks confirm figures exist, are distinct from one another, and are
not drawn over empty categories. They cannot tell you a figure is drawn on the wrong scale — a
mis-scaled colour map passes all three. That is what the `[REVIEW]` items are for, and they are worth
actually looking at.

**Things to check by eye in the `[REVIEW]` figures:**

| Figure | Healthy | Suspicious |
| --- | --- | --- |
| `consensus_T1_vs_T2` | strong diagonal spanning the full category range | a handful of scattered points — the arms are not agreeing |
| `norm_consensus_T1_vs_T2` | bright diagonal on a dark field | a mostly dark matrix — no reproducible categories |
| `state_mu_K_*_arm_*` | a broad, filled cloud roughly ±3 in both axes. **Do not expect separated clusters** — this is the *continuous* within-type state, not the discrete categories | axes far wider than ±3 (an outlier rescaling the plot — cosmetic, see below), or a cloud collapsed to a point |
| `SC_K_*` | most categories above zero, MMIDAS curves near or above the t-type reference | curves hugging zero; or an x-axis spanning a single value, which means the figure is broken rather than the model |
| `classAcc_RF_K_*` | read the **t-types** group — that is MMIDAS vs PCA on reference labels | the `T Categories` groups classify the model's own labels, so ~95% there is near-circular and proves little |
| `conf_*` | tight diagonal with faint off-diagonal detail | large square blocks (reference types collapsing together), or pure black-and-white with no intermediate shades (a plotting-scale bug, not a model result) |
| `state_mu_arm_0_c_*` | highlighted category is a coherent coloured cluster with the traversal path running through it | a single dot, or an invisible highlight |

### Known figure quirks in the distributed run

The figures shipped with this workspace were reviewed by eye. Three things look worse than they are,
and are documented here so you can tell them apart from a genuine problem in your own run.

**1. `state_mu_K_96_arm_0` is unreadable, and this is cosmetic.** A single cell sits at roughly
(9.8, 58) in state space. Matplotlib autoscales to include it, so the y-axis runs to 60 while the
other 22,364 cells occupy about ±3 — they compress into a flat smear along the bottom. The same
figure for arm 1 has no such outlier and shows the expected cloud, and the two arms are otherwise
comparable, which is how you can confirm the smear is a plotting artifact rather than a degenerate
state variable.

To check this on your own run: compare the two arms, and read the axis limits. If one arm's axes are
an order of magnitude wider than the other's, you are looking at an outlier, not a model failure. The
underlying values are in `model.tar.gz` if you want to re-plot with clipped limits.

Worth stressing, because the natural reading is the wrong one: a diffuse cloud here is the **correct**
result. The continuous state captures within-type variation, so it is not supposed to separate into
discrete groups. The discrete structure lives in the categorical variable, and you inspect it in
`consensus_T1_vs_T2` and `norm_consensus_T1_vs_T2` — both of which are clean for this run, showing a
sharp diagonal across all 96 categories with negligible off-diagonal mass.

**2. The `conf_*` heatmaps look mostly white.** Rows are normalised to fractions, so a strong
diagonal leaves everything else pale by construction. That is the intended appearance. What matters
is that intermediate shades are present and the diagonal is continuous; blocks of mid-tone off the
diagonal are real biology — closely related reference t-types being confused with one another. A
`conf_*` figure with *no* intermediate shades at all, only white and full-saturation cells, would
indicate a plotting-scale problem instead.

**3. `classAcc_RF_K_96` shows MMIDAS well below PCA on the t-types group** (about 62% against 91%).
The figure is drawn correctly; this is the run's actual result and the one advisory check that does
not pass. See
[Where the results differ from the published analysis](#where-the-results-differ-from-the-published-analysis).
Note the error bars on the MMIDAS bars are wide (roughly ±5 points), which is the fold-to-fold spread
that makes this metric advisory rather than a fidelity gate.

For reference, the figures that need no such caveat in this run are `consensus_T1_vs_T2_K_96`,
`norm_consensus_T1_vs_T2_K_96`, `state_mu_K_96_arm_1`, and `SC_K_96` — the last showing both arms'
silhouette curves above the t-type PCA baseline across almost the whole range.

### Using it on your own data

On a different dataset there is no published result to compare against, so **treat this notebook as a
starting point for your own validation rather than a pass/fail gate.** Most of it still applies; the
reference comparison does not.

**Keep as-is — none of this depends on the dataset.** All 35 `plumbing` checks: that every expected
file exists, that shapes agree between the `.h5ad`, the manifests and the model, that all three
stages consumed the same inputs (Stage 0 lineage), that manifests agree with one another, and that
figures are distinct rather than redrawn over empty categories. This is the part that tells you the
workflow ran correctly, and it is the same question on any dataset.

**Update to describe your run.** `CONFIG["expected"]` mirrors the inputs you actually passed, so it
has to be re-stated: `n_selected_genes`, `neuronal_classes`, `remove_clusters`, `n_categories`,
`k_select_thr`, `n_selected_cats`, and `kegg_toml_supplied`. The two model-quality floors,
`min_populated_frac` and `min_retained_cpm_frac`, are generic and can stay. Leaving stale ALM/VISp
values here produces failures that say nothing about your run.

**Retire or replace.** `CONFIG["reference"]` and **Stage 6 — Reference comparison** are entirely
Mouse ALM/VISp: `shape`, `n_ttype`, `model_order 92`, `avg_consensus_min`, `pruning_rounds`, and the
three `ttype_acc_*` values all come from the authors' notebooks. Expect the `fidelity` checks to fail
or be meaningless, and disregard the three-way verdict, which is written around them. The
`ttype_acc_*` checks additionally assume curated reference cell-type labels exist to classify
against; if your data has none, they cannot be computed at all.

**What to put in their place.** Substitutes for a published reference, roughly in order of effort:

- **Reproducibility across runs.** Run the same configuration twice and compare `model_order`,
  `avg_consensus` and populated-category counts. Bit-identical results are not expected on a GPU, but
  a stable `model_order` and consensus is the strongest evidence available without a reference.
- **Internal consistency.** `avg_consensus` at or above `k_select_thr`,
  `k_selection_met_threshold` true, and `n_populated_categories` close to `model_order`. These are
  already checked and are meaningful on any dataset.
- **Known biology.** If you have curated labels for even a subset of cells, populate the accuracy
  checks with your own baseline instead of the reference figures. If you have none, the `[REVIEW]`
  figure table above still applies — the healthy-versus-suspicious patterns are properties of the
  model, not of this dataset.

---

## Citation and Credit

MMIDAS is the work of its original authors. If you use these workflows in your research, please cite the original publication:

> Marghi, Y., Gala, R., Baftizadeh, F. & Sümbül, U. Joint inference of discrete cell types and continuous type-specific variability in single-cell datasets with MMIDAS. *Nature Computational Science* **4**, 706–722 (2024). https://doi.org/10.1038/s43588-024-00683-8

**Authors:** Yeganeh Marghi, Rohan Gala, Fahimeh Baftizadeh, and Uygar Sümbül (Allen Institute for Brain Science).

**Original code:** These workflows are built directly on the authors' reference implementation, released by the Allen Institute at **https://github.com/AllenInstitute/MMIDAS**. The `mmidas` Python package and the underlying training, evaluation, clusterability, and state-traversal routines that these WDLs call are the work of the repository's authors and contributors — **Yeganeh Marghi ([@ymarghi](https://github.com/ymarghi))** and **Rohan Gala ([@rhngla](https://github.com/rhngla))**. The five workflow steps here mirror the sequence of the original repository's tutorial notebooks (data preparation → training → evaluation → clusterability → state-traversal analysis).

The scientific method, algorithms, and original software are the intellectual work of these authors and are described in full in the paper and repository above. This workspace only provides WDL wrappers and an example workflow around their published tool; it does not reproduce the publication. Please consult the paper itself for the complete methodology and results, and refer to the [original repository's LICENSE](https://github.com/AllenInstitute/MMIDAS/blob/main/LICENSE) for the terms governing reuse of the MMIDAS software.

## Additional Resources

- **Original publication:** Marghi, Y., Gala, R., Baftizadeh, F. & Sümbül, U. *Joint inference of discrete cell types and continuous type-specific variability in single-cell datasets with MMIDAS.* [Nature Computational Science 4, 706–722 (2024)](https://doi.org/10.1038/s43588-024-00683-8).
- **Original MMIDAS code (Allen Institute):** [github.com/AllenInstitute/MMIDAS](https://github.com/AllenInstitute/MMIDAS) — the reference implementation these workflows are built upon.
- MMIDAS source scripts and Docker image: [warp-tools](https://github.com/broadinstitute/warp-tools) (`3rd-party-tools/mmidas/`).
- WARP repository: [broadinstitute/warp](https://github.com/broadinstitute/warp).
- For Terra-specific documentation and support, see [Terra Support](https://support.terra.bio/hc/en-us).

### Docker image provenance

The `mmidas` image bundles two independently-maintained pieces of code:

1. The **pipeline scripts** (`01_data_prep.py` … `05_state_traversal.py`), which live in
   `warp-tools/3rd-party-tools/mmidas/` and are version-controlled and reviewed there.
2. The **`mmidas` Python package** — a fork of
   [AllenInstitute/MMIDAS](https://github.com/AllenInstitute/MMIDAS) containing the model itself
   (`cpl_mixvae.py`, `eval.py`, `utils/`). Several pipeline behaviours are implemented only here:
   the pruning loop and its `min_con` stop condition, `K_selection`, and the `.h5ad` loader that
   supplies gene identifiers for KEGG pathway mapping.

The second piece is installed from a **pinned revision** of
[jessicaway/MMIDAS](https://github.com/jessicaway/MMIDAS), so any published image can be rebuilt
from source by anyone. `docker_build.sh` resolves the tag in `MMIDAS_GIT_REF` (currently `warp-v2`)
to the immutable commit it points at, fetches that commit's release tarball, and records both:

| Where | What |
| --- | --- |
| Image labels | `MMIDAS_GIT_URL`, `MMIDAS_GIT_REF`, `MMIDAS_SHA` |
| `docker_versions.tsv` | a second column, `<ref>@<sha>`, next to each image tag |

To see exactly which model code an image contains:

```bash
docker inspect --format '{{json .Config.Labels}}' us.gcr.io/broad-gotc-prod/mmidas:<tag>
```

**To pick up new MMIDAS changes:** commit and push them to the fork, move or add a tag, then
rebuild with `./docker_build.sh --mmidas-ref <tag>` (a full 40-character commit SHA also works).
The build fails before doing any work if the ref does not resolve or the tarball is not fetchable.

The image the workflows currently point at is
`us.gcr.io/broad-gotc-prod/mmidas:1.0.0-0.1.0-1787578739`, built at `warp-v2`. Images listed in
`docker_versions.tsv` with the revision `unrecorded` predate this pinning and cannot be rebuilt from
source.

The fork carries four changes on top of upstream `0963ca7`: packaging so the subpackages install, a
headless-plotting fix and a threshold-comparison fix in `K_selection`, the gene-identifier fix in
`load_data` that KEGG pathway mapping depends on, and added per-epoch entropy and per-round consensus
logging (stdout only). The `load_data` fix is a candidate for upstreaming to AllenInstitute/MMIDAS.

---

## Contact Information

- For workspace questions and feedback, email the Broad pipelines team at warp-pipelines-help@broadinstitute.org.
- You can also contact the Terra team from the Terra main menu. When submitting a request, it is helpful to include your Project ID, workspace name, Bucket ID, Submission ID, and Workflow ID, and any relevant log information.

---

## License

**Copyright Broad Institute, 2026 | BSD-3**
All code provided in this workspace is released under the WDL open source code license (BSD-3) (full license text at https://github.com/broadinstitute/warp/blob/develop/LICENSE). Note however that the programs called by the scripts may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running these tools.

---

## Workspace Change Log

| Date | Change |
| --- | --- |
| 2026-09-03 | Reviewed the distributed run's figures by eye and documented three that look worse than they are, most importantly that a diffuse `state_mu` cloud is the correct result rather than a sign of collapse. |
| 2026-09-02 | Documented run-to-run variation with two independent runs of the same configuration (`model_order` 96 and 89, `avg_consensus` 0.9075 and 0.9692, bracketing the published 92 and 0.939). Added a section on MMIDAS's multimodal (transcriptomic + electrophysiology) capability, which these workflows do not implement, and what porting it would involve. |
| 2026-08-28 | Clarified that `MMIDAS_output_validation.ipynb` validates against the authors' published Mouse ALM/VISp results, and documented which parts transfer to a different dataset and what to use in place of the reference comparison. |
| 2026-08-27 | `MMIDAS_Train` `1.3.0`: `n_epoch_p` default is now **1,000**, the configuration this workspace's example outputs were produced with (~13.5 h, ~$13). The published analysis used 10,000 (~4.5 days, ~$110); set it explicitly to reproduce that. Corrected the inputs table, which had understated `max_prun_it` and described `min_con` as halting pruning. Time and cost figures are now measured and stated once. |
| 2026-08-26 | Documented dataset-specific tuning (`tau` scaling with `n_categories`, `max_prun_it` bounding reachable `model_order`, runtime/cost scaling, no-resume exposure). Recorded the validated result against the published analysis and where it differs. Added the validation-notebook section. |
| 2026-08-21 | Figure fixes in `MMIDAS_Analyze` (`1.1.2`): confusion matrices row-normalised before plotting, accuracy bar chart widened, state-traversal category labels made unique, highlight palette no longer collides with the background. |
| 2026-08-11 | Corrected training defaults (`tau`, `x_drop`, `max_prun_it`) so the published `model_order` is reachable and the categorical posterior commits at `n_categories = 120`. Pinned the Docker build to a tagged MMIDAS revision so images are reproducible. |
| 2026-08-06 | Fixed KEGG pathway mapping, checkpoint selection in `03a_evaluate.py`, state-traversal category selection, and the silhouette figure. Added review fields to `evaluation_results.json`. |
| 2026-06-24 | Initial MMIDAS example workspace documentation. |
