# Mitochondria Pipelines

This directory contains pipelines for processing mitochondrial DNA (mtDNA) data, including variant calling, coverage analysis, and merging coverage metrics based on the mtSwirl v2.5_MongoSwirl_Single pipeline and mtSwirl scripts for running in Terra. All original WDL and scripts can be found the [mtSwirl GitHub repository](https://github.com/rahulg603/mtSwirl).

This pipeline was originally released as part of the manuscript: Nuclear genetic control of mitochondrial DNA copy number and heteroplasmy in humans, which can be found at [Nature](https://www.nature.com/articles/s41586-023-06426-5). 

If you use these resources in your work, please cite as Gupta et al. 2023 Nature:
Gupta, R., Kanai, M., Durham, T.J. et al. Nuclear genetic control of mtDNA copy number and heteroplasmy in humans. Nature, in press. https://doi.org/10.1038/s41586-023-06426-5.

The pipelines are designed to handle hg38-aligned input files and produce high-quality outputs for downstream analysis.

---

## **Mitochondria Pipeline**

### Overview

The Mitochondria Pipeline processes mtDNA data from whole-genome sequencing (WGS) BAM or CRAM files. It generates high-quality variant calls, coverage metrics, and other mitochondrial-specific outputs. The pipeline includes steps for alignment, variant calling, and coverage analysis.

### Key Features

- **Input Compatibility**: Supports hg38-aligned BAM or CRAM files.
- **Variant Calling**: Produces VCF files with SNP and INDEL calls specific to mtDNA.
- **Coverage Analysis**: Generates base-level coverage metrics for mtDNA and NUMTs.
- **Self-Reference Generation**: Creates self-references for improved variant calling accuracy.
- **NUMT Analysis**: Optionally computes NUMT coverage metrics for quality control.

### Inputs

- **WGS BAM/CRAM**: Aligned to hg38.
- **Reference Files**:
    - `ref_fasta`, `ref_fasta_index`, `ref_dict`: Reference genome files.
    - `mt_fasta`, `mt_fasta_index`, `mt_dict`: Mitochondrial reference files.
    - `mt_interval_list`, `nuc_interval_list`: Interval lists for mtDNA and NUMTs.
- **Blacklisted Sites**: Files to exclude problematic regions.
- **Optional Inputs**:
    - `max_read_length`: Read length for optimization.
    - `vaf_filter_threshold`: Threshold for filtering low VAF sites.

### Outputs

- **Variant Files**:
    - `final_vcf`: Final VCF of mtDNA SNPs and INDELs.
    - `self_ref_vcf`: VCF generated using self-references.
- **Coverage Metrics**:
    - `final_base_level_coverage_metrics`: Base-level coverage metrics for mtDNA.
    - `numt_base_level_coverage` (optional): NUMT coverage metrics.
- **Intermediate Files**:
    - `subset_bam`, `subset_bai`: BAM/BAI files for mtDNA.
    - `self_reference_fasta`: Generated self-reference FASTA file.

---

## **mt_coverage_merge Pipeline**

### Overview

The `mt_coverage_merge` workflow takes per-sample mtDNA VCF files and coverage metrics produced by the Mitochondria Pipeline and merges them into a single annotated cohort-wide callset. The pipeline is designed to scale to hundreds of thousands of samples.

The workflow proceeds in the following stages:
1. **(Optional) Sample subsetting** — filter the input data table to a specific sample list
2. **TSV preprocessing** — join auxiliary metadata (ancestry, age, WGS coverage) into a master sample sheet
3. **Coverage DB construction** — build an HDF5 coverage database from per-sample coverage files
4. **VCF ingestion (sharded + merged)** — import per-sample VCFs in parallel shards and merge via a multi-round fan-in tree
5. **Finalization (sharded)** — impute homoplasmic-reference genotypes from the coverage DB and apply artifact filters, run in parallel across sample shards
6. **Annotation** — add cohort-wide variant statistics, QC filters, and population/haplogroup allele frequencies; run twice in parallel (all samples and QC-filtered samples)

### Inputs

| Input | Type | Description |
|---|---|---|
| `full_data_tsv` | File | Per-sample data table with paths to VCF and coverage files from the Mitochondria Pipeline |
| `sample_list_tsv` | File? | Optional list of sample IDs to subset to |
| `coverage_tsv` | File | Supplementary QC metrics (mean WGS coverage, contamination, collection date) |
| `ancestry_tsv` | File | Predicted genetic ancestry per sample |
| `dob_tsv` | File | Date of birth per sample |
| `wgs_median_coverage_tsv` | File | WGS median coverage per sample |
| `step3_output_bucket` | String | GCS bucket for intermediate MT tarballs |
| `annotated_output_bucket` | String | GCS bucket for final annotated outputs |
| `vcf_col_name` | String | Column name in the data table containing the VCF path (default: `final_vcf`) |
| `combined_mt_name` | String | Name for the combined output MT (default: `combined_vcf`) |

Key sharding parameters (with defaults):

| Parameter | Default | Description |
|---|---|---|
| `step3_shard_size` | 2,500 | Samples per VCF ingest shard |
| `step3_merge_fanin` | 10 | Merge fan-in per round |
| `step3_shard_n_partitions` | 192 | Hail partitions per shard MT |
| `finalize_shard_size` | 25,000 | Samples per finalize shard |
| `finalize_shard_n_partitions` | 256 | Hail partitions per finalize shard MT |
| `finalize_union_n_partitions` | 1,000 | Hail partitions for the final unioned MT |

### Outputs

| Output | Type | Description |
|---|---|---|
| `processed_tsv` | File | Master sample sheet after metadata joins |
| `output_coverage_db` | File | `coverage_db.tar.gz` — HDF5 coverage database + summary TSV |
| `combined_vcf` | File | `tar.gz` of the finalized, unioned cohort MatrixTable |
| `annotated_output_gcs_path` | String | GCS path to annotated outputs (all samples) |
| `filt_annotated_output_gcs_path` | String | GCS path to annotated outputs (QC-filtered samples only) |

---

### Tasks

#### `subset_data_table` *(optional)*
Filters `full_data_tsv` to only the rows matching sample IDs in `sample_list_tsv`. If no sample list is provided, the full table is passed through unchanged.

#### `process_tsv_files`
Joins the sample data table with `coverage_tsv`, `ancestry_tsv`, `dob_tsv`, and `wgs_median_coverage_tsv` on sample ID. Renames and derives columns: `ancestry_pred` → `pop`, collection date − date of birth → `age`. The output `processed_tsv` is used as the master sample sheet for all downstream tasks.

#### `annotate_coverage`
Builds an HDF5 coverage database (`coverage.h5`) from the per-sample coverage files listed in `processed_tsv`. Per-position summary statistics (mean, median, fraction over 100× and 1000×) are computed using numpy. Outputs `coverage_db.tar.gz` containing `coverage.h5` and `coverage_summary.tsv`.
- **Docker**: `aou-mito-coverage-db:1.0.0`

#### `make_vcf_shards_from_tsv`
Partitions the sample list into TSV shards of `step3_shard_size` samples each, producing one TSV file per shard.

#### `build_vcf_shard_mt` *(scattered)*
Imports one shard of per-sample VCF files into a Hail MatrixTable, writes it as a `.tar.gz` to `step3_output_bucket`, and returns the GCS path.
- **Docker**: `aou-mito-hail-processing:1.0.0`

#### `make_mt_merge_groups` + `merge_mt_shards` *(3 rounds, scattered)*
Merges shard MTs in a multi-round fan-in tree (`step3_merge_fanin = 10`). Each round groups the current MT tarballs and merges each group via `multi_way_union_mts`. Three rounds reduce ~200+ shards to a single merged MT. Intermediate tarballs are stored in `step3_output_bucket`.
- **Docker**: `aou-mito-hail-processing:1.0.0`

#### `shard_mt_by_samples`
Splits the fully merged MT into per-sample-shard MT tarballs of `finalize_shard_size` samples each, stored in `step3_output_bucket`.
- **Docker**: `aou-mito-hail-processing:1.0.0`

#### `finalize_mt_with_covdb` *(scattered)*
For each sample shard: reads the shard MT and the HDF5 coverage DB, imputes homoplasmic-reference genotypes (missing entries where coverage exceeds threshold are set to hom-ref), and applies the artifact-prone site filter. Processes positions in blocks to avoid repeated full-matrix scans. Outputs a finalized shard MT tarball.
- **Docker**: `aou-mito-hail-processing:1.0.0`

#### `union_mt_shards`
Unions all finalized sample-shard MTs into a single cohort-wide MT via `union_cols`, writes the result as a tarball to `step3_output_bucket`.
- **Docker**: `aou-mito-hail-processing:1.0.0`

#### `add_annotations` *(called twice: `annotated` and `filt_annotated`)*
Runs `add_annotations.py` on the final unioned MT to add cohort-wide variant statistics, per-population and per-haplogroup allele frequencies, QC filters, and pathogenicity annotations. Called twice in parallel:
- `annotated` (`keep_all_samples = true`): all samples retained
- `filt_annotated` (`keep_all_samples = false`): samples failing QC (contamination, copy number, etc.) removed

Outputs are written directly to timestamped subdirectories under `annotated_output_bucket`.
- **Docker**: `aou-mito-hail-processing:1.0.0`
