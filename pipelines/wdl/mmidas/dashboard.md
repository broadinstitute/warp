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

| Input | Type | Default | Description |
| --- | --- | --- | --- |
| `preprocessed_h5ad` | File | — | Output of `MMIDAS_DataPrep`, or your own contract-compliant `.h5ad` |
| `run_augmenter` | Boolean | `false` | Whether to train the optional data augmenter first |
| `n_categories` | Int | `120` | Upper-bound number of categories before pruning |
| `n_arm` | Int | `2` | Number of coupled encoder arms |
| `state_dim` | Int | `2` | Continuous within-type state dimension |
| `latent_dim` | Int | `10` | Low-dimensional embedding dimension |
| `training_mode` | String | `MSE` | Reconstruction loss (`MSE` or `ZINB`) |
| `n_epoch` / `n_epoch_p` | Int | `10000` / `1000` | Epochs before pruning / per pruning round |
| `max_prun_it` | Int | `14` | Maximum pruning iterations |
| `min_con` | Float | `0.99` | Inter-arm consensus at which pruning stops early. Pruning continues while any surviving category is below this, up to `max_prun_it` rounds |
| `k_select_thr` | Float | `0.95` | Consensus threshold used to recommend `model_order` |
| `batch_size` | Int | `5000` | Mini-batch size |
| `seed` | Int | `0` | Random seed |
| `train_gpu` | Int | `0` | Set to `1` to attach a GPU (see note below) |

> **GPU note.** Training is much faster on a GPU — the measured benchmarks above assume one. Set `train_gpu = 1` and that is all: the `TrainMixVAE` task's runtime block already declares `gpuCount` and `gpuType: "nvidia-tesla-t4"` and passes `--cuda` to the training script. There is **no GPU setting to enable in Terra** for workflow submissions — Terra's Cloud Environment has GPU options, but those apply to interactive notebooks and RStudio, not to Cromwell workflow tasks.
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
model reconstructs entirely through the continuous state variable, ignoring the discrete categories
it is supposed to be learning. We hit exactly this: with `tau = 0.1` at `n_categories = 120` (a
value that is roughly correct for `n_categories = 15`), a full training run produced

- `avg_consensus` **0.026** instead of ~0.97,
- **11** populated categories out of 120,
- a `model_order` that was meaningless.

The same run with `tau = 0.005` gave `avg_consensus` 0.969 and 89 of 89 categories populated. The
failure is silent unless you look: nothing errors, and `model_order` still comes back a plausible
number.

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

Total epochs are `n_epoch + max_prun_it × n_epoch_p`, and at ~0.9 s/epoch on an `nvidia-tesla-t4`
(~$1.00/hr) that dominates everything else in the pipeline:

| `n_epoch_p` | Total epochs | Time | Cost |
| --- | --- | --- | --- |
| 1,000 (`MMIDAS_Train.staged_validation.json`) | 52,000 | ~14 h | **~$14** |
| 10,000 (`MMIDAS_Train.json`, reference value) | 430,000 | ~5 d | **~$120** |

**Run the staged configuration first.** On the example data the cheap run reached
`avg_consensus 0.969` and `model_order 89`, matching the published analysis — the 10× longer
configuration was not needed. Round-by-round consensus in a full-length run plateaued by round 2
and then moved only within noise for 17 more rounds.

### 4. There is no resume — protect against losing a long run

`TrainMixVAE` writes `model.tar.gz` only when the task **completes**. If it is aborted (a cost cap,
a timeout, a manual cancel), the intermediate checkpoints are lost with the VM and the run must
start over. `preemptible` is already `0` so the VM will not be reclaimed mid-run, but:

- **Set Terra cost caps above the expected spend** (≥ ~$25 for the staged config, ≥ ~$150 for the
  full-length one). We lost ~2 days and ~$50 of a full-length run to a forgotten cap.
- 4–5 days is within the usual 7-day GCP task ceiling, but only just. Confirm your project does not
  impose a shorter limit before launching the full-length configuration.

### 5. Smaller things worth knowing

| Parameter / behaviour | What to know |
| --- | --- |
| `min_con` | **Reporting only.** The consensus-based pruning stop is commented out upstream, so pruning always runs the full `max_prun_it`. Do not expect `min_con` to halt anything. |
| `k_select_thr` | If no checkpoint reaches it, `K_selection` returns nothing and `Evaluate` falls back to the un-pruned checkpoint. Check `k_selection_met_threshold` in `evaluation_results.json` — `false` means `model_order` came from a fallback, not a selection. |
| `kegg_toml` | Optional. Omit it and `n_pathways` is 0 with no pathway figures — expected, not a failure. Supply it and confirm `n_pathways > 0`; zero pathways *with* a `kegg_toml` means gene-name lookup failed. |
| `htree_file` | Optional; enables taxonomy ordering in stage 03c. |
| `n_selected_cats` | Capped at the number of *populated* categories, so the manifest may report fewer than you asked for. |
| Run-to-run variation | Training is **not** bit-reproducible even with a fixed `seed`, because GPU reductions are non-deterministic. Two runs of the identical configuration diverged enough that the two arms swapped which one converged better. Expect `model_order` to move by a few categories between runs. |
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

Measured on Terra with the example dataset (22,365 cells x 5,032 genes). Training dominates; the
other two stages are minor by comparison.

### MMIDAS_DataPrep

| Dataset | Cells | Genes | Time | Cost |
| --- | --- | --- | --- | --- |
| Mouse ALM-VISp (example) | 22,365 | 5,032 | ~25 min | < $1 |

Note this stage needs `mem_size = 48` GiB: it holds the full 22,439 x 45,768 count matrix in memory
before subsetting to the selected genes. The earlier default of 32 GiB was not enough.

### MMIDAS_Train

Measured at ~0.9 s/epoch on an `nvidia-tesla-t4` at ~$1.00/hr, with total epochs
`n_epoch + max_prun_it x n_epoch_p`.

| Config | GPU | n_epoch | n_epoch_p | max_prun_it | Total epochs | Time | Cost |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `MMIDAS_Train.staged_validation.json` | T4 | 10,000 | 1,000 | 42 | 52,000 | ~14 h | **~$14** |
| `MMIDAS_Train.json` (reference values) | T4 | 10,000 | 10,000 | 42 | 430,000 | ~5 d | **~$120** |

Run the staged configuration first — see
[Validate cheaply before paying for the full run](#validate-cheaply-before-paying-for-the-full-run).
CPU-only training was not benchmarked; it is impractical at these epoch counts.

### MMIDAS_Analyze

| Dataset | Time | Cost |
| --- | --- | --- |
| Mouse ALM-VISp (example) | ~1 h | < $2 |

CPU only, no GPU. Cheap enough to re-run freely, which matters because it is the stage you re-run
when figures or downstream analysis change without retraining.

For more information about controlling Cloud costs, see [this article](https://support.terra.bio/hc/en-us/articles/360029748111).

---

## Fidelity to the original MMIDAS analysis

These workflows are a **port**, not a reimplementation or an improvement: the goal is that running
them in Terra does what running the authors' notebooks does. The
[MMIDAS repo](https://github.com/AllenInstitute/MMIDAS) ships its notebooks with outputs saved, so
the authors' results for the example dataset are recorded and can be compared against directly.

**Reference values for Mouse ALM/VISp:**

| Quantity | Reference | Source |
| --- | --- | --- |
| matrix shape | 22,365 × 5,032 | `1_data_prep.ipynb` |
| reference t-types | 115 | `2_train.ipynb` data summary |
| pruning rounds | 42 | `3_evaluation.ipynb` (checkpoints `after_pruning_1..42`) |
| `model_order` | **92** | `3_evaluation.ipynb`; hardcoded in notebooks 4 and 5 |
| `avg_consensus` | 0.939 (test cells) / 0.954 (K-selection) | `3_evaluation.ipynb` |

**Measured result for the validated run in this workspace:**

| Quantity | Reference | This workspace | |
| --- | --- | --- | --- |
| matrix shape | 22,365 × 5,032 | 22,365 × 5,032 | match |
| reference t-types | 115 | 115 | match |
| pruning rounds | 42 | 42 | match |
| `model_order` | 92 | **89** | within tolerance (3 pruning rounds) |
| `avg_consensus` | 0.939 / 0.954 | **0.9692** | above the published value |
| populated categories | 92 (all) | **89 of 89** (both arms) | all populated |
| K-selection met `k_select_thr` | yes | **yes** | no fallback |

`MMIDAS_output_validation.ipynb` checks a run against these under **Stage 6 — Reference
comparison**; see [Validating a run](#validating-a-run--mmidas_output_validationipynb) for how to
point it at your own submissions and what the check labels mean. Where the results *do* differ, see
[Where the results differ from the published analysis](#where-the-results-differ-from-the-published-analysis).

### Where the training defaults come from

The reference notebooks pass only `n_categories`, `state_dim`, `n_arm` and `latent_dim` to
`cpl_mixVAE.init_model()` and let everything else fall to that function's defaults. `MMIDAS_Train`'s
defaults follow those, with training-loop values (`n_epoch`, `n_epoch_p`, `max_prun_it`) taken from
`tutorials/train_mixvae.py`, since `2_train.ipynb`'s `n_epoch = 10` is a walkthrough demo rather than
the configuration that produced the published model.

Two of these matter enough to call out:

- **`max_prun_it` bounds which answers are reachable.** Pruning removes one category per round, so
  the smallest `model_order` a run can produce is `n_categories - max_prun_it`. Reaching the
  reference `model_order` of 92 from `n_categories = 120` requires 28 rounds; the default is 42, the
  number the reference used. Lowering it below 28 puts the published answer *outside the search
  space*, and no amount of training will recover it.
- **`tau` is tied to `n_categories`.** `init_model` documents it as "usually equals to
  1/n_categories" — about 0.0083 at `n_categories = 120`. The default here is the library's 0.005.
  Do not carry over `tutorials/train_mixvae.py`'s `0.1`: that script's `n_categories` default is 15,
  where 0.1 is roughly 1/K. At K=120 it leaves the categorical posterior far too soft to reach the
  consensus the reference reports.

`min_con` is **reporting only**. The consensus-based pruning stop is commented out in
`mmidas/cpl_mixvae.py::train` upstream, so pruning always runs the full `max_prun_it`; every
reference invocation behaves this way. Do not enable it — with `min_con = 0.99` it halts pruning as
soon as one category is reproducible, which stops short of the depth the published results rely on.

### Validate cheaply before paying for the full run

`max_prun_it = 42` at `n_epoch_p = 10000` is roughly **4.5 days on an `nvidia-tesla-t4` (~$110)**,
measured at 0.90 s/epoch and $1.00/hr. `TrainMixVAE` runs with `preemptible: 0` and has no resume
path, so a failure late in that run is expensive.

`example_inputs/MMIDAS_Train.staged_validation.json` is the same configuration with
`n_epoch_p = 1000` — **13 hours, about $13**. It differs from the reference config in that one field
only. Use it first:

| Config | Epochs | Time | Cost |
| --- | --- | --- | --- |
| `MMIDAS_Train.staged_validation.json` | 52,000 | 13 h | ~$13 |
| `MMIDAS_Train.json` (reference) | 430,000 | 4.5 days | ~$110 |

The staged run still has `max_prun_it = 42`, so `model_order = 92` is inside the search space
(the floor is `120 - 42 = 78`). What it buys is a cheap answer to the one open question — whether the
categorical posterior commits once `tau` is scaled to `n_categories`.

**There is an even earlier decision point.** The categorical posterior's behaviour is visible in the
`TrainMixVAE` log from the first epochs, before any pruning happens:

```
Entropy: -9.0249 (uniform=-9.5750)
```

`uniform` is the value that column takes when both arms' categorical posteriors are flat, i.e. the
discrete latent carries no information. If `Entropy` is still pinned near `uniform` at the end of the
pre-pruning phase (the first `n_epoch = 10000` epochs, ~2.5 h, ~$2.50), the temperature change did
not work and there is no point continuing — kill the run rather than paying for 42 pruning rounds on
a collapsed model. If it has moved substantially toward 0, let it finish.

### Known divergences from the notebooks

Two differences are unavoidable in a generic workflow and are deliberate:

| Divergence | Why |
| --- | --- |
| **State-traversal category selection.** `5_state_traversal.ipynb` hardcodes `selected_c = [80, 119, 1, 13, 92, 25, 55, 110, 62, 69, 31]`. The workflow instead selects the `n_selected_cats` most-populated categories. | A workflow cannot reproduce a hand-picked list chosen by inspection. Selecting by population at least guarantees the figures show categories the model actually uses. |
| **Train/test split seeding.** `2_train.ipynb` calls `get_loaders` without a `seed`, so its split is random and unrecoverable. The workflow seeds it (`seed = 0`). | An unseeded split makes a workflow non-reproducible run to run. This is why an exact numerical match to the reference is not expected, and why the validation notebook compares `model_order` within a tolerance and `avg_consensus` against a range. |

One further note: the pipeline scripts import `load_data` from `mmidas.utils.data_tools`, whereas the
reference notebooks import it from `mmidas.utils.dataloader`. Both read `adata.X` and the `var` index
identically, so the expression matrix handed to training is the same; `dataloader` additionally
derives `cluster_id` / `c_onehot` helpers that the notebooks use for their own plots and the
workflow scripts compute where needed.

### Where the results differ from the published analysis

The validated run reproduces the reference within tolerance, but three things do not match exactly.
None indicates a broken port; all three are worth knowing before you present results.

**1. `model_order` 89 vs the published 92.** Three pruning rounds apart. Two causes, neither
fixable: the reference train/test split was unseeded, and GPU training is not bit-reproducible even
at a fixed seed. Two runs of our own identical configuration diverged enough that the two encoder
arms swapped which one converged better. Expect a few categories of movement between runs, and do
not treat any single `model_order` as *the* answer.

**2. t-type classification accuracy is further below the PCA baseline than in the reference.**
This is the one metric that does not match well:

| | PCA-100 | MMIDAS-10 | gap |
| --- | --- | --- | --- |
| Reference (`4_clusterability.ipynb`) | ~84.5% | ~73.5% | **~0.11** |
| This run | 90.9% | 66.2% | **0.247** |

Our PCA baseline does *better* than the reference's and our MMIDAS embedding does *worse*, so the
gap is roughly double. Every comparison points the same direction as the reference — PCA wins on the
115 reference t-types, the 10-dimensional MMIDAS embedding wins decisively on MMIDAS's own
categories (94.7% vs 64.8%) — so this reads as a magnitude difference in a downstream metric rather
than a fidelity failure. Two caveats on the comparison itself: this metric has ~0.05 fold-to-fold
spread, and the reference numbers were **read off a bar chart by eye** because
`4_clusterability.ipynb` prints no values. The validation notebook reports this as *advisory* for
exactly those reasons.

Note also that a gap here is expected by design: a 10-dimensional embedding is not attempting to
beat a 100-component PCA basis at recovering 115 reference labels.

**3. The accuracy bar chart has fewer groups than the reference's.** The reference plots three label
sets — t-types, *Merged t-types*, and MMIDAS T Categories. The workflow plots t-types and the
per-arm MMIDAS categories, omitting the merged-t-type group (reference t-types collapsed down to
`model_order` groups). The two groups that do appear match the reference's layout and direction.

---

## Validating a run — `MMIDAS_output_validation.ipynb`

The workspace includes a notebook that checks a completed set of runs end to end. Point it at your
three submissions, run it top to bottom, and it reports what passed, what failed, and which figures
still need a human eye. It is the fastest way to know whether a run is trustworthy before you build
analysis on top of it.

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

| Kind | Question | On failure |
| --- | --- | --- |
| `plumbing` | Did the workflow execute correctly? Files present, shapes consistent, manifests mutually agreeing, all stages consuming the same inputs. | The port is broken. Fix before interpreting anything. |
| `fidelity` | Does this run reproduce the published analysis? Compared against values recorded in the reference notebooks. | The port runs but not the way the authors ran it. |
| `advisory` | Is the model any good? Soft metrics with wide run-to-run spread, or comparisons against figures read by eye. | Informational. Never fatal. |

A run that reproduces the reference is a success even if advisory items complain — and a run with
better-looking numbers that does *not* reproduce the reference is not.

**Result for the validated example run:** 47/48 — `plumbing 35/35`, `fidelity 8/8`, `advisory 4/5`,
verdict *"the workflow executed correctly and matches the reference analysis within tolerance"*. The
single advisory failure is the t-type accuracy gap described above.

**What it cannot tell you.** The checks confirm figures exist, are distinct from one another, and are
not drawn over empty categories. They cannot tell you a figure is drawn on the wrong scale — a
mis-scaled colour map passes all three. That is what the `[REVIEW]` items are for, and they are worth
actually looking at: an earlier round of this workspace shipped confusion matrices that rendered as
black-and-white noise (raw counts plotted against a 0–1 colour scale) and every automated check
passed.

**Things to check by eye in the `[REVIEW]` figures:**

| Figure | Healthy | Suspicious |
| --- | --- | --- |
| `consensus_T1_vs_T2` | strong diagonal spanning the full category range | a handful of scattered points — the arms are not agreeing |
| `norm_consensus_T1_vs_T2` | bright diagonal on a dark field | a mostly dark matrix — no reproducible categories |
| `state_mu_K_*_arm_*` | visibly separated clusters | one undifferentiated blob, or far fewer groups than `model_order` |
| `SC_K_*` | most categories above zero, MMIDAS curves near or above the t-type reference | curves hugging zero; or an x-axis spanning a single value, which means the figure is broken rather than the model |
| `classAcc_RF_K_*` | read the **t-types** group — that is MMIDAS vs PCA on reference labels | the `T Categories` groups classify the model's own labels, so ~95% there is near-circular and proves little |
| `conf_*` | tight diagonal with faint off-diagonal detail | large square blocks (reference types collapsing together), or pure black-and-white with no intermediate shades (a plotting-scale bug, not a model result) |
| `state_mu_arm_0_c_*` | highlighted category is a coherent coloured cluster with the traversal path running through it | a single dot, or an invisible highlight |

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
from source by anyone. `docker_build.sh` resolves the tag in `MMIDAS_GIT_REF` (default `warp-v1`)
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

Images built before `warp-v1` are listed in `docker_versions.tsv` with the revision `unrecorded`.
Those were built by copying a local, uncommitted working tree and **cannot be reproduced**; do not
treat them as a known quantity.

The fork carries four changes on top of upstream `0963ca7`: packaging so the subpackages install, a
headless-plotting fix and a threshold-comparison fix in `K_selection`, the gene-identifier fix in
`load_data`, and the restored consensus-based pruning stop in `train`. The last two are candidates
for upstreaming to AllenInstitute/MMIDAS.

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
| 2026-08-26 | Documented dataset-specific tuning (`tau` scaling with `n_categories`, `max_prun_it` bounding reachable `model_order`, runtime/cost scaling, no-resume exposure). Recorded the validated result against the published analysis and where it differs. Added the validation-notebook section. Replaced placeholder cost estimates with measured Terra numbers. |
| 2026-08-21 | Figure fixes in `MMIDAS_Analyze` (`1.1.2`): confusion matrices row-normalised before plotting, accuracy bar chart widened, state-traversal category labels made unique, highlight palette no longer collides with the background. |
| 2026-08-11 | Corrected training defaults to match the reference analysis (`tau`, `x_drop`, `n_epoch_p`, `max_prun_it`); reverted the `min_con` pruning stop to upstream behaviour. Pinned the Docker build to a tagged MMIDAS revision so images are reproducible. |
| 2026-08-06 | Fixed KEGG pathway mapping, checkpoint selection in `03a_evaluate.py`, state-traversal category selection, and the silhouette figure. Added review fields to `evaluation_results.json`. |
| 2026-06-24 | Initial MMIDAS example workspace documentation. |
