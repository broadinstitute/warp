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

The `mt_coverage_merge` pipeline processes mitochondrial coverage data by combining multiple input files, annotating the coverage, and generating a final output. It is designed to handle large-scale datasets and produce merged coverage metrics for downstream analysis.

### Key Features

- **Input Merging**: Combines multiple coverage-related TSV files into a single processed file.
- **Annotation**: Annotates coverage data with additional metadata (e.g., ancestry, date of birth).
- **Hail Integration**: Uses Hail for efficient processing and annotation of large datasets.

### Inputs

- **Coverage Data**:
    - `coverage_tsv`: A TSV file containing coverage metrics.
    - `ancestry_tsv`: A TSV file with ancestry predictions.
    - `dob_tsv`: A TSV file with date of birth data.
- **Sample Table Fields**:
    - Arrays of strings representing various genomic metrics and files (e.g., `n_liftover_r2_spanning_complex`, `self_mt_aligned_bai`, etc.).
- **Docker Image**:
    - `hail_docker_img`: Specifies the Docker image for running Hail-based tasks.

### Outputs

- **Processed TSV**:
    - `processed_tsv`: A merged and annotated TSV file containing coverage metrics and metadata.
- **Hail Table**:
    - `output_ht`: A compressed Hail Table for downstream analysis.

### Workflow Steps

1. **Generate TSV File**:
    - Combines input arrays into a single TSV file.
2. **Process TSV Files**:
    - Merges the combined TSV with additional metadata (e.g., ancestry, date of birth).
    - Calculates participant age and other derived metrics.
3. **Annotate Coverage**:
    - Uses Hail to annotate the processed TSV file and generate a Hail Table.

### Runtime Requirements

- **Docker Images**:
    - `hail_docker_img`: For Hail-based tasks.
- **Resources**:
    - Memory and CPU requirements vary by task, with high-memory configurations for large-scale processing.

---

## Usage

### Inputs

Provide the required input files and parameters in a JSON or YAML format. Example:

```json
{
  "coverage_tsv": "path/to/coverage.tsv",
  "ancestry_tsv": "path/to/ancestry.tsv",
  "dob_tsv": "path/to/dob.tsv",
  "hail_docker_img": "hailgenetics/hail:0.2.67"
}