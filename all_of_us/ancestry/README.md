# Ancestry Pipelines
The following pipelines are used to calculate ancestry:

## determine_hq_sites_intersection.wdl

### Background
Identifies a set of high-quality (HQ) variant sites that are common between a pre-defined training set and a new dataset, and then filter both the training set and the new data to only include genotypes at these shared HQ sites. This is often a preparatory step for downstream analyses like ancestry inference or population structure modeling.

### Overview
#### Inputs

* training_vcf_bgz: Full training VCF with sample genotypes.

* training_vcf_so_bgz: Sites-only version of the training VCF (used for speed).

* ordered_vcf_shards_in: List of input data VCF shards to be processed.

#### Step 1 – Filter Sites in Input Data Shards
- For each input VCF shard:
  - Use `GATK SelectVariants` to filter:
    - Only biallelic SNPs
    - Minor allele frequency (AF) > 0.001
    - No more than 1% missing genotypes
    - Exclude filtered variants
    - Optionally restrict to provided genomic intervals (e.g., exome regions)
  - Output: a sites-only VCF per shard.

#### Step 2 – Intersect Filtered Sites with Training HQ Sites
- For each filtered VCF shard:
  - Intersect the sites-only VCF with the training set's HQ sites-only VCF.
  - Output: an intersected sites-only VCF per shard.

#### Step 3 – Merge Intersected Sites Across Shards
- Combine all intersected VCFs into a single, merged sites-only VCF.
- This file contains the set of HQ sites that are present in both:
  - The training set
  - The input dataset

#### Step 4 – Filter Training VCF to Intersection Sites
- Filter the full training VCF to include only the sites present in the merged intersection VCF.
- Output: a filtered training VCF with genotypes at shared HQ sites.

#### Step 5 – Filter Input Data VCF Shards to Intersection Sites
- For each original input VCF shard:
  - Filter it to only retain variants at the merged HQ intersection sites.
- Output: a filtered VCF for each shard.

#### Step 6 – Merge Filtered Input Data Shards
- Combine all filtered input VCF shards into one final merged VCF.
- Output: a single, complete dataset VCF containing only HQ sites shared with the training data.

## run_ancestry
This workflow trains a PCA-based ancestry classifier using a reference dataset (e.g., HGDP) and applies the model to infer ancestry for new samples.

---

#### Inputs

| Name | Description |
|------|-------------|
| `hq_variants_intersection` | Sites-only VCF of HQ sites shared between training and input data |
| `hq_variants_intersection_idx` | Index for `hq_variants_intersection` |
| `merged_vcf_shards` | Merged VCF of input dataset at HQ sites |
| `merged_vcf_shards_idx` | Index for `merged_vcf_shards` |
| `filtered_training_set` | Training dataset VCF at HQ sites |
| `filtered_training_set_idx` | Index for `filtered_training_set` |
| `hgdp_metadata_file_in` *(optional)* | Sample metadata for the training set (e.g. HGDP) |
| `final_output_prefix` | Prefix for naming all outputs |
| `other_cutoff_in` *(optional)* | Probability cutoff for labeling as "Other" |
| `num_pcs` | Number of principal components to use (default: 16) |

---

#### Step 1 – `create_hw_pca_training`
- Perform PCA on the `filtered_training_set`
- Join sample metadata to PCA scores
- Output: PCA scores, loadings, population labels, eigenvalues

#### Step 2 – `call_ancestry`
- Project input samples (`merged_vcf_shards`) into the PCA space using training loadings
- Apply Random Forest classifier trained on PCA features of training data
- Classify ancestry based on prediction probabilities (with optional “Other” cutoff)
- Output: Predicted ancestry labels, probabilities, classifier model, Hail tables

#### Step 3 – `plot_ancestry`
- Generate interactive HTML plots of ancestry predictions using the first two PCs
- Two versions: one with original predictions, one with "Other" classification applied

---

#### Outputs

| Output | Description |
|--------|-------------|
| `results_tsv` | TSV of ancestry predictions, probabilities, PCA features |
| `results_ht` | Hail Table of projected PCA scores with sample IDs |
| `results_loadings_ht` | Hail Table of PCA loadings with allele frequencies |
| `classifier_pkl` | Pickled Random Forest classifier |
| `pred_plot` | PCA plot with original ancestry predictions |
| `pred_oth_plot` | PCA plot with "Other" predictions based on cutoff |
| `training_pca_labels_ht_tsv` | PCA scores + population labels for training samples |
| `eigenvalues_txt` | Eigenvalues used in PCA |