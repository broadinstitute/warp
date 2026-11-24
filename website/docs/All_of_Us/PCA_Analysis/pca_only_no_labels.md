---
sidebar_position: 1
slug: /All_of_Us/PCA_Analysis/pca_only_no_labels
title: PCA (No Labels)
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.0.0](https://github.com/broadinstitute/warp/blob/develop/all_of_us/PCA/pca_only_no_labels.changelog.md) | August, 2026 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

<!-- TODO: Add a pca_only_no_labels workflow diagram here once available, e.g.:
![pca_only_no_labels workflow diagram](./pca_only_no_labels_diagram.png)
-->

## Introduction to the pca_only_no_labels workflow

[`pca_only_no_labels`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/PCA/pca_only_no_labels.wdl) is a WDL workflow that performs Hardy-Weinberg equilibrium–normalized Principal Component Analysis (PCA) on genomic variant data using [Hail](https://hail.is/). It is designed for exploratory analysis of population structure **without** requiring pre-existing population labels, and it produces both tabular results (per-sample scores and eigenvalues) and a set of visualization plots.

The workflow accepts one or more BGZF-compressed VCF files (with Tabix indices), concatenating them when more than one is provided. It can optionally subset the cohort to a single predicted ancestry before training. It then runs Hail's `hwe_normalized_pca` for a configurable number of components and emits a scores TSV, an eigenvalues TSV, a scree plot, and per-PC-pair scatter, hexbin, and interactive 3D density plots. Because this workflow does not use population labels, all samples receive a placeholder "No label" designation in the plots.

## Quickstart table

| Pipeline Feature | Description | Source |
| :--: | :-- | :--: |
| Analysis type | Hardy-Weinberg normalized PCA (unlabeled) | [Hail `hwe_normalized_pca`](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.hwe_normalized_pca) |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Genomic reference sequence | GRCh38 | |
| Data input file format | BGZF-compressed VCF + Tabix index (`.tbi`) | |
| Data output file format | TSV, PNG, interactive HTML | |
| Primary software | Hail, bcftools, matplotlib, plotly | [Hail](https://hail.is/), [bcftools](https://samtools.github.io/bcftools/), [matplotlib](https://matplotlib.org/), [Plotly](https://plotly.com/python/) |

## Set-up

### pca_only_no_labels installation and requirements

The workflow code can be downloaded by cloning the [WARP GitHub repository](https://github.com/broadinstitute/warp). For the latest release, please see the [pca_only_no_labels changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/PCA/pca_only_no_labels.changelog.md).

The pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant workflow management system, such as in a [Terra](https://app.terra.bio) environment.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `hq_sites_vcf_files` | One or more BGZF-compressed VCF files containing variant data. When more than one is provided, they are assumed to be in chromosomal order and are concatenated. | Array[File] |
| `hq_sites_vcf_indices` | Tabix index (`.tbi`) files corresponding to `hq_sites_vcf_files`. Must be non-empty and the same length as `hq_sites_vcf_files`. | Array[File] |
| `final_output_prefix` | Prefix applied to all output filenames. | String |
| `num_pcs` | Number of principal components to compute. | Int |
| `min_vcf_partitions_in` | *(Optional)* Minimum number of partitions for VCF import. Default: `100`. | Int? |
| `alpha` | *(Optional)* Scatter-plot point opacity. Default: `0.18`. | Float |
| `pc_pairs` | *(Optional)* List of PC pairs to plot. If omitted, defaults to PC1 vs PC2 and PC3 vs PC4. If provided, **only** the given pairs are plotted (the list replaces, rather than extends, the defaults). JSON form: `[{"left":5,"right":6},{"left":7,"right":8}]`. | Array[Pair[Int, Int]]? |
| `ancestry_list` | *(Optional)* TSV (with header) of ancestry predictions covering all samples. Must contain a `research_id` column (sample ID) and an `ancestry_pred_other` column (ancestry label); other columns are ignored. **Provide together with `ancestry`.** | File? |
| `ancestry` | *(Optional)* Which `ancestry_pred_other` value to subset to (e.g. `"eur"`). **Provide together with `ancestry_list`.** | String? |

:::note
`ancestry_list` and `ancestry` must be supplied **together**. Providing exactly one (but not both) is an error; providing neither runs PCA on the full cohort.
:::

## pca_only_no_labels tasks and tools

The [pca_only_no_labels workflow](https://github.com/broadinstitute/warp/blob/develop/all_of_us/PCA/pca_only_no_labels.wdl) calls a series of tasks to concatenate inputs (when needed), train the PCA model, compute per-PC variance proportions, and generate plots.

Overall, the workflow:

1. Validates the input VCF/index arrays.
2. Concatenates per-chromosome VCFs (only when more than one is provided).
3. Trains the HWE-normalized PCA model and exports scores and eigenvalues.
4. Computes each PC's proportion of variance among the computed PCs.
5. Generates a scree plot.
6. Generates scatter, hexbin, and interactive 3D density plots for each requested PC pair.

To see specific tool parameters, select the task WDL link in the table; then view the `command {}` section of the task in the WDL script. To view or use the exact tool software, see the task's Docker image which is specified in the task WDL `runtime` section as `docker:`.

| Task name and WDL link | Tool | Software | Description |
| --- | --- | --- | --- |
| [ConcatenateChromosomalVcfs](https://github.com/broadinstitute/warp/blob/develop/all_of_us/PCA/pca_only_no_labels.wdl) | concat, index | [bcftools](https://samtools.github.io/bcftools/) | *(Conditional — only when more than one VCF is provided.)* Concatenates the per-chromosome BGZF VCFs into a single file and creates a `.tbi` index. |
| [create_hw_pca_training](https://github.com/broadinstitute/warp/blob/develop/all_of_us/PCA/pca_only_no_labels.wdl) | `hwe_normalized_pca` | [Hail](https://hail.is/) | Imports the VCF into Hail, optionally subsets columns to a single ancestry, runs HWE-normalized PCA (loadings not computed) for `num_pcs` components, and exports per-sample scores and eigenvalues as TSVs. |
| [compute_pct_variance](https://github.com/broadinstitute/warp/blob/develop/all_of_us/PCA/pca_only_no_labels.wdl) | Python (pandas) | [pandas](https://pandas.pydata.org/) | Computes, for each PC, `100 × eigenvalue / sum(computed eigenvalues)`. |
| [plot_scree](https://github.com/broadinstitute/warp/blob/develop/all_of_us/PCA/pca_only_no_labels.wdl) | matplotlib | [matplotlib](https://matplotlib.org/) | Plots each PC's variance proportion vs PC index as a scree plot. |
| [plot_pca](https://github.com/broadinstitute/warp/blob/develop/all_of_us/PCA/pca_only_no_labels.wdl) | matplotlib, plotly | [matplotlib](https://matplotlib.org/), [Plotly](https://plotly.com/python/) | *(Scattered over PC pairs.)* Produces a scatter plot, a log-scale hexbin density plot, and an interactive 3D density HTML surface for each requested PC pair. |

### 1. Input validation

The workflow fails early (via `Utilities.ErrorWithMessage`) if `hq_sites_vcf_files` and `hq_sites_vcf_indices` are empty or differ in length.

### 2. ConcatenateChromosomalVcfs (conditional)

Runs **only when more than one** VCF file is provided; a single input VCF is used directly. The per-chromosome BGZF VCFs are concatenated with `bcftools concat` and indexed with `bcftools index`. The output basename defaults to `<final_output_prefix>_autosomes.vcf.gz`.

### 3. create_hw_pca_training

Imports the (concatenated or single) BGZF VCF into Hail with the requested minimum partitions. When `ancestry_list` and `ancestry` are supplied, the MatrixTable columns are subset to the samples whose `ancestry_pred_other` equals the requested ancestry **before** any training; sample IDs are matched exactly between the table's `research_id` and the VCF sample name (`s`). This subsetting hard-fails with a descriptive error if the ancestry matches no rows, if any listed sample is absent from the MatrixTable, or if there is zero overlap. The task then runs Hail's `hwe_normalized_pca` (loadings not computed) for `num_pcs` components and exports per-sample scores (columns `s`, `PC1`…`PC{num_pcs}`) and the eigenvalues as TSVs.

### 4. compute_pct_variance

Reads the eigenvalues and computes, for each PC, `100 × eigenvalue / sum(computed eigenvalues)`.

:::note Interpretation
For each retained PC, the proportion of variance is calculated as that PC's eigenvalue divided by the sum of eigenvalues among the `k` computed PCs (`k` = `num_pcs`). Each value is therefore that PC's share of the **retained eigenvalue mass** — it sums to 100% over the selected PCs by construction — and is **not** the true proportion of total genome-wide variance, since only the top `k` components were computed (it does not normalize by the trace / all eigenvalues). This same metric is used consistently for the plot axis labels and the scree plot. It is labeled "Proportion of variance among computed PCs (%)" throughout, with column header `Variance_Proportion_Among_Computed_PCs_Pct`.
:::

### 5. plot_scree

Reads the `compute_pct_variance` output directly (single source of truth) and plots the per-PC value vs PC index as a scree plot. The "elbow" where the curve levels off indicates how many PCs capture meaningful structure.

### 6. plot_pca (scattered over PC pairs)

Runs once per plotted PC pair and produces, for each pair, a scatter plot, a hexbin density plot on a logarithmic color scale (viridis colormap), and an interactive 3D density surface as self-contained HTML. Axis labels include the per-PC variance proportion (e.g. `PC1 (12.34%)`). Each requested pair is validated in-task: the two PCs must differ, and both must exist among the computed PCs (otherwise the task fails with a clear message).

## Outputs

| Output variable name | Filename, if applicable | Output format and description |
| --- | --- | --- |
| `training_pca_labels_ht_tsv` | `<final_output_prefix>_training_pca.tsv` | TSV of PCA scores for all samples (columns `s`, `PC1`…`PC{num_pcs}`). |
| `training_pca_eigenvalues_tsv` | `<final_output_prefix>_training_pca_eigenvalues.tsv` | TSV of PCA eigenvalues. |
| `training_pca_scree_plot` | `<final_output_prefix>_scree.png` | Scree plot PNG (proportion of variance among computed PCs, %). |
| `training_pca_scatter_plots` | `<final_output_prefix>_<pc1>_<pc2>_scatter.png` | One scatter-plot PNG per plotted PC pair, ordered to match `pc_pairs` (or the defaults). |
| `training_pca_hexbin_plots` | `<final_output_prefix>_<pc1>_<pc2>_hexbin.png` | One log-scale hexbin density PNG per plotted PC pair. |
| `training_pca_3d_density_interactive_plots` | `<final_output_prefix>_<pc1>_<pc2>_3d_density.html` | One interactive 3D density HTML surface per plotted PC pair. |

## Runtime requirements

Runtime attributes below reflect the defaults defined in the WDL. `ConcatenateChromosomalVcfs` runs only when more than one VCF is provided.

| Task | Docker image | Memory | CPU | Disk |
| --- | --- | :--: | :--: | --- |
| ConcatenateChromosomalVcfs | `mgibio/bcftools-cwl:1.12` | 128 GB | 16 | 1.5 TB HDD |
| create_hw_pca_training | `hailgenetics/hail:0.2.134-py3.11` | 512 GB | 48 | 2 TB SSD |
| compute_pct_variance | `us.gcr.io/broad-gotc-prod/warp-tools:2.6.1` | 16 GB | 2 | 250 GB HDD |
| plot_scree | `faizanbashir/python-datascience:3.6` | 8 GB | 2 | 100 GB HDD |
| plot_pca | `faizanbashir/python-datascience:3.6` | 16 GB | 2 | 500 GB HDD |

For `create_hw_pca_training`, the task reserves ~50 GB for the OS/Python process and assigns the remainder to the Spark driver heap via `SPARK_DRIVER_MEMORY`.

## Important notes

- **Ancestry subsetting** requires `ancestry_list` and `ancestry` to be supplied together, with `ancestry` matching an `ancestry_pred_other` value exactly (case-sensitive). Because runs are per-ancestry, set `final_output_prefix` per run (e.g. include the ancestry) so outputs don't collide, and ensure the `research_id` values match the VCF sample names exactly — mismatches hard-fail with counts and example IDs rather than silently corrupting results.
- **PC pairs:** plotting defaults to PC1 vs PC2 and PC3 vs PC4; set `pc_pairs` to plot a different set (it replaces the defaults). Both PCs in a pair must be within `1..num_pcs`.
- **Network dependency:** `plot_pca` installs `plotly==5.18.0` via `pip` at task runtime (it is not baked into the plotting image) to render the interactive 3D HTML. The task will fail if the compute environment has no PyPI egress.
- **No population labels:** all samples receive a placeholder "No label" designation, since this workflow does not use population labels.

## Versioning

All pca_only_no_labels pipeline releases are documented in the [pca_only_no_labels changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/PCA/pca_only_no_labels.changelog.md). The All of Us versioning convention (for example, `aou_9.0.0`) reflects the release version in which the pipeline was used.

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
