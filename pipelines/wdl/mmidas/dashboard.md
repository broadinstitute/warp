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
```

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

`avg_consensus` should also be at or above `k_select_thr`. A run with high `model_order` and near-zero `avg_consensus` has not found reproducible categories; it has failed to train its discrete latent. In the training log, watch the per-epoch `Entropy` against the `uniform=` value printed next to it — an `Entropy` that stays pinned at `uniform` while the reconstruction loss falls means the categorical variable never committed and no amount of pruning will fix it.

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

> **GPU note.** Training is much faster on a GPU. Set `train_gpu = 1` to attach an `nvidia-tesla-t4`. On Terra, also enable GPU in the runtime/quota settings for your project.

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

<!-- TODO: Populate with measured benchmarks from Terra runs. -->

Benchmarks below are placeholders to be filled in once measured on Terra. Training time depends heavily on GPU vs. CPU, `n_epoch`, `max_prun_it`, and dataset size.

### MMIDAS_DataPrep

| Dataset | Cells | Genes | Time | Cost $ |
| --- | --- | --- | --- | --- |
| Mouse ALM-VISp (example) | _TBD_ | ~1,252 | _TBD_ | _TBD_ |

### MMIDAS_Train

| Dataset | GPU | n_epoch | max_prun_it | Time | Cost $ |
| --- | --- | --- | --- | --- | --- |
| Mouse ALM-VISp (example) | T4 | 10000 | 14 | _TBD_ | _TBD_ |
| Mouse ALM-VISp (example) | CPU only | 10000 | 14 | _TBD_ | _TBD_ |

### MMIDAS_Analyze

| Dataset | Time | Cost $ |
| --- | --- | --- |
| Mouse ALM-VISp (example) | _TBD_ | _TBD_ |

For more information about controlling Cloud costs, see [this article](https://support.terra.bio/hc/en-us/articles/360029748111).

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

| Date | Change | Author |
| --- | --- | --- |
| _TBD_ | Initial MMIDAS example workspace documentation. | _TBD_ |
