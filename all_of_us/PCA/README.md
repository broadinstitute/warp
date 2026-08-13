# PCA Analysis Pipeline

The following pipeline performs Principal Component Analysis (PCA) on genomic variant data without population labels.

## pca_only_no_labels

### Background

This WDL workflow performs Hardy-Weinberg equilibrium normalized PCA on genomic variant data from BGZ-compressed VCF files, using [Hail](https://hail.is/) for large-scale processing. It is designed for exploratory analysis of population structure without requiring pre-existing population labels, and it produces both tabular results (scores and eigenvalues) and a set of visualization plots.

Key characteristics:
- Accepts one or more BGZ-compressed VCF files (with indices); concatenates them when more than one is provided.
- Optionally subsets the cohort to a single predicted ancestry before training.
- Performs Hardy-Weinberg equilibrium normalized PCA for a configurable number of components.
- Emits a scores TSV, an eigenvalues TSV, a scree plot, and per-PC-pair scatter, hexbin, and interactive 3D density plots.

### Inputs

- `Array[File] hq_sites_vcf_files` – BGZ-compressed VCF files containing variant data (assumed to be in chromosomal order when more than one).
- `Array[File] hq_sites_vcf_indices` – Index files (`.tbi`) corresponding to `hq_sites_vcf_files`. Must be the same length as `hq_sites_vcf_files` and non-empty.
- `String final_output_prefix` – Prefix for all output filenames.
- `Int num_pcs` – Number of principal components to compute.
- `Int? min_vcf_partitions_in` – Optional minimum number of partitions for VCF import (default: 100).
- `Float alpha` – Scatter-plot point opacity (default: 0.18).
- `Array[Pair[Int, Int]]? pc_pairs` – Optional list of PC pairs to plot. If omitted, defaults to PC1 vs PC2 and PC3 vs PC4. If provided, **only** the given pairs are plotted (include every pair you want — it replaces, not extends, the defaults). JSON form: `[{"left":5,"right":6},{"left":7,"right":8}]`.
- `File? ancestry_list` – Optional TSV (with header) of ancestry predictions covering all samples. Must contain a `research_id` column (sample ID) and an `ancestry_pred_other` column (ancestry label); other columns are ignored. **Provide together with `ancestry`.**
- `String? ancestry` – Which `ancestry_pred_other` value to subset to (e.g. `"eur"`). **Provide together with `ancestry_list`.**

Providing exactly one of `ancestry_list` / `ancestry` (but not both) is an error. Providing neither runs PCA on the full cohort.

### Workflow steps

#### Input validation
Fails early (via `Utilities.ErrorWithMessage`) if `hq_sites_vcf_files` and `hq_sites_vcf_indices` are empty or differ in length.

#### Step 1. ConcatenateChromosomalVcfs (conditional)
Runs **only when more than one** VCF file is provided; a single input VCF is used directly.
- Concatenates the per-chromosome BGZ VCFs into a single file with `bcftools concat` and creates a `.tbi` index with `bcftools index`.
- Output basename defaults to `<final_output_prefix>_autosomes.vcf.gz`.

#### Step 2. create_hw_pca_training
- Imports the (concatenated or single) BGZ VCF into Hail with the requested minimum partitions.
- **Optional ancestry subsetting:** when `ancestry_list` + `ancestry` are supplied, the MatrixTable columns are subset to the samples whose `ancestry_pred_other` equals the requested ancestry, **before** any training. Sample IDs are matched exactly between the table's `research_id` and the VCF sample name (`s`). This step hard-fails with a descriptive error if the ancestry matches no rows, if any listed sample is absent from the MatrixTable, or if there is zero overlap.
- Runs Hail's `hwe_normalized_pca` (loadings not computed) for `num_pcs` components.
- Exports per-sample scores (columns `s`, `PC1`…`PC{num_pcs}`) and the eigenvalues as TSVs.

#### Step 3. compute_pct_variance
- Reads the eigenvalues and computes, for each PC, `100 * eigenvalue / sum(computed eigenvalues)`.
- **Note on interpretation:** this normalizes by the sum of the *computed* (top-`num_pcs`) eigenvalues, not the total variance (trace / all eigenvalues). Each value is therefore a PC's share of the *retained* eigenvalue mass — it sums to 100% over the selected PCs by construction and is **not** the true proportion of total variance explained. It is labeled "Proportion of variance among computed PCs (%)" throughout, with column header `Variance_Proportion_Among_Computed_PCs_Pct`.

#### Step 4. plot_scree
- Reads the `compute_pct_variance` output directly (single source of truth) and plots the per-PC value vs PC index as a scree plot. The "elbow" where the curve levels off indicates how many PCs capture meaningful structure.

#### Step 5. plot_pca (scattered over PC pairs)
Runs once per plotted PC pair and produces, for each pair:
- A **scatter** plot (`<prefix>_<pc1>_<pc2>_scatter.png`).
- A **hexbin** density plot on a logarithmic color scale using the viridis colormap (`<prefix>_<pc1>_<pc2>_hexbin.png`).
- An **interactive 3D density** surface as self-contained HTML (`<prefix>_<pc1>_<pc2>_3d_density.html`).

Axis labels include the per-PC variance proportion (e.g. `PC1 (12.34%)`). Each requested pair is validated in-task: `pc1` and `pc2` must differ, and both must exist among the computed PCs (otherwise the task fails with a clear message).

### Outputs

- `File training_pca_labels_ht_tsv` – TSV of PCA scores for all samples (columns `s`, `PC1`…`PC{num_pcs}`).
- `File training_pca_eigenvalues_tsv` – TSV of PCA eigenvalues.
- `File training_pca_scree_plot` – Scree plot PNG (proportion of variance among computed PCs, %).
- `Array[File] training_pca_scatter_plots` – One scatter-plot PNG per plotted PC pair, ordered to match `pc_pairs` (or the defaults). Files are named by pair, e.g. `<prefix>_1_2_scatter.png`.
- `Array[File] training_pca_hexbin_plots` – One log-scale hexbin density PNG per plotted PC pair.
- `Array[File] training_pca_3d_density_interactive_plots` – One interactive 3D density HTML per plotted PC pair.

### Runtime requirements

**ConcatenateChromosomalVcfs** (only when >1 VCF):
- Docker: `mgibio/bcftools-cwl:1.12`
- Memory: 128 GB · CPU: 16 · Disk: 1.5 TB HDD

**create_hw_pca_training:**
- Docker: `hailgenetics/hail:0.2.134-py3.11`
- Memory: 512 GB · CPU: 48 · Disk: 2 TB SSD
- Reserves ~50 GB for the OS/Python and assigns the remainder to the Spark driver heap via `SPARK_DRIVER_MEMORY`.

**compute_pct_variance:**
- Docker: `us.gcr.io/broad-gotc-prod/warp-tools:2.6.1`
- Memory: 16 GB · CPU: 2 · Disk: 250 GB HDD

**plot_scree:**
- Docker: `faizanbashir/python-datascience:3.6`
- Memory: 8 GB · CPU: 2 · Disk: 100 GB HDD

**plot_pca:**
- Docker: `faizanbashir/python-datascience:3.6`
- Memory: 16 GB · CPU: 2 · Disk: 500 GB HDD

### Usage notes

- **Ancestry subsetting** requires `ancestry_list` and `ancestry` to be supplied together, with `ancestry` matching an `ancestry_pred_other` value exactly (case-sensitive). Because runs are per-ancestry, set `final_output_prefix` per run (e.g. include the ancestry) so outputs don't collide, and ensure the `research_id` values match the VCF sample names exactly — mismatches hard-fail with counts and example IDs rather than silently corrupting results.
- **PC pairs:** plotting defaults to PC1 vs PC2 and PC3 vs PC4; set `pc_pairs` to plot a different set (it replaces the defaults). Both PCs in a pair must be within `1..num_pcs`.
- **Network dependency:** `plot_pca` installs `plotly==5.18.0` via `pip` at task runtime (it is not in the plotting image) to render the interactive 3D HTML. The task will fail if the compute environment has no PyPI egress.
- All samples receive a placeholder "No label" designation since this workflow does not use population labels.
