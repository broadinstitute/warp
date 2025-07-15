# PCA Analysis Pipeline
The following pipeline is used to perform Principal Component Analysis (PCA) on genomic variant data without population labels.

## pca_only_no_labels
#### Background

This WDL workflow performs Principal Component Analysis on genomic variant data from BGZ-compressed VCF files.
It uses Hail for efficient processing of large-scale genomic datasets and generates PCA scores along with visualization plots.
The workflow is designed for exploratory analysis of population structure in genomic data without requiring pre-existing population labels.

Key characteristics:
- Processes BGZ-compressed VCF files with indices
- Performs Hardy-Weinberg equilibrium normalized PCA
- Generates configurable number of principal components
- Creates scatter plots for visualization of PC1 vs PC2
- Handles large datasets with configurable partitioning
- Outputs both tabular results and visualization plots

#### Inputs
Analysis Parameters:
- `File full_bgz` – BGZ-compressed VCF file containing variant data
- `File full_bgz_index` – Index file for the BGZ-compressed VCF
- `String final_output_prefix` – Prefix for all output filenames
- `Int num_pcs` – Number of principal components to compute
- `Int? min_vcf_partitions_in` – Optional minimum number of partitions for VCF processing (default: 200)

#### Step 1. create_hw_pca_training
- Input Validation: Ensures BGZ file and index are properly formatted and accessible
- VCF Import: Loads the BGZ-compressed VCF file into Hail with specified minimum partitions
- PCA Computation: Performs Hardy-Weinberg equilibrium normalized PCA on genotype data
- Data Export: Flattens PCA scores and exports results as a TSV file
- Resource Management: Uses high-memory configuration for processing large datasets

Technical Details:
- Uses Hail's `hwe_normalized_pca` function for population structure analysis
- Processes genotype (GT) field from VCF data
- Exports eigenvalues and scores but not loadings for efficiency
- Configured with 123 GB memory and 16 CPUs for large-scale processing

#### Step 2. plot_pca
- Data Loading: Reads the PCA results TSV file into a pandas DataFrame
- Score Parsing: Converts string representations of score arrays back to numeric arrays
- Label Assignment: Assigns "No label" to all samples since no population labels are available
- Visualization: Creates scatter plots comparing specified principal components
- Plot Generation: Saves publication-ready PNG plots with proper labeling and legends

Technical Details:
- Validates PC indices to ensure they are positive (1-indexed)
- Uses matplotlib for high-quality scientific plotting
- Applies rainbow color scheme for potential future categorical data
- Generates scatter plots with 2-point markers for clarity

#### Outputs

- `File training_pca_labels_ht_tsv` – TSV file containing PCA scores for all samples
- `File training_pca_labels_tsv_plots` – PNG plot file showing PC1 vs PC2 scatter plot

#### Runtime Requirements

**create_hw_pca_training task:**
- Docker: `hailgenetics/hail:0.2.67`
- Memory: 123 GB
- CPU: 16 cores
- Disk: 500 GB HDD

**plot_pca task:**
- Docker: `faizanbashir/python-datascience:3.6`
- Memory: 16 GB
- CPU: 2 cores
- Disk: 500 GB HDD

#### Usage Notes

- The workflow is optimized for large-scale genomic datasets
- BGZ files are processed with localization_optional flag for improved performance
- Default partitioning of 200 can be adjusted based on dataset size
- Plot generation defaults to PC1 vs PC2 but can be configured for other component pairs
- All samples receive "No label" designation since this workflow doesn't use population labels