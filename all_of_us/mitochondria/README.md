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

The workflow takes in a TSV file containing sample names, as well as output metrics and output files files from the Mitochondria Pipeline workflow. Its goal is to produce a final, comprehensively annotated mitochondrial variant callset. It does this by:
- Preparing a master sample information file.
- Merging all individual coverage files into a single coverage dataset.
- Merging all individual VCF files, using the coverage data to improve accuracy.
- Running an extensive annotation and QC process on the final merged VCFs.

A key feature of this pipeline is that the final annotation step is run twice in parallel: once keeping all samples (`annotated`) and once filtering out low-quality samples (`filt_annotated`), producing two distinct final outputs.

### Task: `subset_data_table`
- **Purpose**: This is an optional first step that filters the input sample TSV (containing output metrics from the Mitochondria Pipeline) to only include a specific list of desired samples.
- **Inputs**:
- `full_data_tsv`: The input TSV containing information for all samples, including paths to their VCF and coverage files produced by the Mitochondria Pipeline
  - `sample_list_tsv` (Optional): A simple text file with one column of sample IDs to keep.

- **Transformation**:
  - A simple Python script using pandas reads both files.
  - It filters the `full_data_tsv` to keep only the rows whose sample ID is present in the `sample_list_tsv`.
  - If `sample_list_tsv` is not provided, it simply passes the `full_data_tsv` through unchanged.

- **Output**:
  - `subset_tsv`: A new TSV file containing data for only the selected samples.

### Task: `process_tsv_files`
- **Purpose**: To perform data wrangling by merging several different information sources into a single, clean, master sample sheet that the downstream Hail scripts can use.

- **Inputs**:
  - `input_tsv`: The (potentially subsetted) data table from the previous step.
  - `coverage_tsv`: A supplementary file containing QC metrics like mean sequencing coverage.
  - `ancestry_tsv`: A file with the predicted genetic ancestry for each sample.
  - `dob_tsv`: A file with the date of birth for each sample.

- **Transformation**:
  - This task runs a pandas script to join these four files together based on the sample ID.
  - It renames columns to match the expected format for the Hail scripts (e.g., `ancestry_pred` becomes `pop`).
  - It calculates a new `age` column by subtracting the date of birth from the biosample collection date.

- **Output**:
  - `processed_tsv`: A single, master TSV file with all the necessary columns (`s`, `coverage`, `final_vcf`, `pop`, `age`, etc.) for the main pipeline.

### Task: `annotate_coverage`
- **Purpose**: To merge all the individual per-base coverage files into a single, unified Hail MatrixTable.
- **Inputs**:
  - `input_tsv`: The master processed_tsv from the previous step, which contains the paths to each sample's coverage file.

- **Transformation**:
  - This task calls the `annotate_coverage.py` script. 
  - This script performs the scalable, **hierarchical merge** of all the individual coverage files listed in the input TSV. 
  - After merging, it calculates summary statistics across all samples, such as the mean and median coverage at each base of the mitochondrial genome. 
  - Finally, it archives the resulting Hail MatrixTable into a compressed tarball (`.tar.gz`).

- **Output**:
  - `output_ht`: A `tar.gz` archive containing the final combined coverage MatrixTable.

### Task: `combine_vcfs`
- **Purpose**: To merge all individual sample VCF files into a single, combined MatrixTable, using the coverage data for improved accuracy.
- **Inputs**:
  - `input_tsv`: The master `processed_tsv`, which contains the paths to each sample's VCF file.
  - `coverage_mt_tar`: The compressed coverage MatrixTable produced by the `annotate_coverage` task.

- **Transformation**:
  - This task calls the `combine_vcfs.py` script.
  - It first unzips the coverage MatrixTable so it can be read by Hail.
  - The Python script then performs the **hierarchical merge** of all the individual VCF files.
  - Critically, it uses the unzipped coverage data to accurately distinguish between sites with no variation (homoplasmic reference) and sites with low sequencing depth (missing data).
  - The resulting merged VCF MatrixTable is archived into a `.tar.gz` file.

- **Output**:
  - `results_tar`: A `tar.gz` archive containing the combined VCF MatrixTable.

### Task: `add_annotations` (called as `annotated` and `filt_annotated`)
- **Purpose**: To perform the final, extensive annotation and quality control of the combined VCF callset.
- **Inputs**:
  - `vcf_mt`: The compressed VCF MatrixTable from the `combine_vcfs` task.
  - `coverage_mt`: The compressed coverage MatrixTable from the `annotate_coverage` task.
  - `coverage_tsv`: The master `processed_tsv`, which is used here to provide sample-level QC stats.
  - `keep_all_samples`: A boolean flag that is the key difference between the two calls.

- **Transformation**:
  - This task calls the comprehensive `add_annotations.py` script. It first unzips the input VCF and coverage MatrixTables.
  - The script then adds dozens of annotations, including:
    - Allele frequencies (overall, per-population, per-haplogroup). 
    - Pathogenicity predictions for tRNA variants. 
    - Quality control filters (`indel_stack`, `common_low_heteroplasmy`, etc.). 
  - The transformation depends on the `keep_all_samples` flag:
    - `annotated` call (`keep_all_samples = true`): The script calculates all QC metrics but does not remove any samples. 
    - `filt_annotated` call (`keep_all_samples = false`): The script removes samples that fail QC filters for contamination, copy number, etc.

- **Output**:
  - `annotated_output_tar` / `filt_annotated_output_tar`: Two separate `.tar.gz` archives, each containing a complete set of output files (MatrixTables, VCFs, text files) for the respective analysis—one with all samples included, and one with only the QC-passed samples.
