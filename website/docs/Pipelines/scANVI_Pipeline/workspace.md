---
sidebar_position: 2
slug: /Pipelines/scANVI_Pipeline/workspace
---

# scANVI Terra Workspace

| Workspace | Maintained by | Questions or Feedback |
| :---: | :---: | :---: |
| [Multiome SCVI and SCANVI](https://app.terra.bio/#workspaces/warp-pipelines/Multiome%20SCVI%20and%20SCANVI) | WARP Pipelines | Please [file an issue in WARP](https://github.com/broadinstitute/warp/issues) |

## Overview

The **Multiome SCVI and SCANVI** Terra workspace is a ready-to-run home for the [scANVI pipeline](./README.md). It comes preloaded with the scANVI WDL workflow configuration and example Multiome 10k PBMC data, so you can run production-scale **cell type label transfer** with no additional setup. The workspace also includes an interactive notebook that walks through the same label-transfer logic step by step, using [scVI](https://docs.scvi-tools.org/en/1.4.1/user_guide/models/scvi.html) and [scANVI](https://docs.scvi-tools.org/en/1.4.1/user_guide/models/scanvi.html) (single-cell ANnotation using Variational Inference).

The workspace offers two entry points:

- **The scANVI WDL workflow** (recommended for production and at-scale runs) — automated, reproducible label transfer via Cromwell/Terra. Start here.
- **The "Transfer Labels Using SCANVI Model" notebook** (for interactive exploration) — the same logic run cell by cell; see [The interactive notebook](#the-interactive-notebook) below.

## Running the scANVI WDL workflow

The workspace ships with a preloaded WDL workflow configuration for running scANVI at scale. The workflow is fully automated and suitable for production use, and mirrors the notebook logic.

1. Open the [Multiome SCVI and SCANVI](https://app.terra.bio/#workspaces/warp-pipelines/Multiome%20SCVI%20and%20SCANVI) workspace on Terra.
2. Navigate to the **Workflows** tab and select the preloaded **scANVI** workflow.
3. Set the inputs. At minimum provide `input_id` plus either the direct file inputs (`gex_h5ad` and `ref_h5ad`, and optionally `atac_h5ad`) or an `input_bucket`. To run **inference only** from a previously trained model, supply `scanvi_model` and set `gpu_count = 0` for a cheaper CPU-only run. See the [Inputs section of the pipeline overview](./README.md) for the full input list and reference requirements.
4. Launch the analysis. Training runs request GPUs by default (`gpu_count = 2`); CPU-only prediction runs set `gpu_count = 0`.

For full documentation — all inputs, outputs, task descriptions, and runtime configuration — see the [scANVI pipeline overview](./README.md).

### Example run costs

Representative Terra costs for scANVI runs across three references (mouse hippocampus, human neocortex, and PBMC) in both training and CPU-inference modes. Costs are from single runs on GCP `us-central1` (training on **2× NVIDIA T4** GPUs; inference CPU-only via `gpu_count = 0`) and will vary with region, machine type, preemption, and dataset. "Cells labelled" is the number of query cells annotated (after PreprocessFilter; for multiome, after the GEX↔ATAC shared-barcode intersection); "Reference" is the annotated `ref_h5ad`.

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
- **Per-cell cost rises as the query shrinks**, because the fixed reference-load / training cost is amortized over fewer labelled cells.
- **Multiome and GEX-only cost about the same** ($3.64 vs $3.38) — the ATAC branch adds little.

## The interactive notebook

The workspace also provides the **"Transfer Labels Using SCANVI Model"** notebook for interactive, step-by-step exploration of the label-transfer process. It focuses specifically on **annotating an ATAC query using a gene expression matrix that already carries cell type labels**, using human Multiome datasets (10k PBMC) processed through the WARP [PeakCalling](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/peak_calling) pipeline alongside the PBMC reference from the [scArches tutorial](https://docs.scarches.org/).

The notebook walks through:

1. **Downloading input files** — retrieves the GEX, ATAC, and reference h5ad files needed for the analysis.
2. **Environment setup** — installs and configures the required Python packages, including [SnapATAC2](https://kzhang.org/SnapATAC2/) and [scvi-tools](https://docs.scvi-tools.org/).
3. **Data filtering and processing** — filters cells, aligns barcodes between the GEX and ATAC modalities, and converts the ATAC cell-by-bin matrix to a gene activity matrix.
4. **Label transfer** — trains SCVI and SCANVI models and transfers cell type labels from the annotated GEX reference onto the unannotated ATAC query, following the [SnapATAC2 label transfer tutorial](https://kzhang.org/SnapATAC2/).

### Running the notebook

The workspace is preloaded with the Multiome 10k PBMC data, so no additional setup is required before launching the notebook.

1. Open the [Multiome SCVI and SCANVI](https://app.terra.bio/#workspaces/warp-pipelines/Multiome%20SCVI%20and%20SCANVI) workspace on Terra.
2. Navigate to the **Analysis** tab.
3. Select the **"Transfer Labels Using SCANVI Model"** notebook.
4. Choose a Virtual Machine and **Enable GPUs** — the SCANVI model training step requires a GPU.
5. Run the notebook cells in sequence to perform the analysis.

:::note GPU requirement
Model training in this notebook requires a GPU-enabled cloud environment. When prompted to select a runtime, make sure **Enable GPUs** is checked before starting the VM.
:::

:::tip Relationship to the WDL pipeline
The notebook walks through the same label-transfer logic implemented in the scANVI WDL pipeline, but interactively — making it a useful companion for understanding the pipeline or for exploratory analysis before running production-scale jobs. See the [scANVI pipeline documentation](./README.md) for running at scale via Cromwell or Terra workflows.
:::
