# HLA Pipelines
The following pipelines are used to generate and aggregate HLA genotyping results from sequencing data.

## HLAGenotyping
#### Background

This WDL workflow performs HLA genotyping using three separate tools—HLA-HD, Polysolver, and OptiType—and generates a consensus genotype call.

Key characteristics:
- Uses GATK to extract and prepare HLA-specific reads.
- Runs HLA-HD, and conditionally runs Polysolver and OptiType when HLA-HD fails to emit a three-field allele.
- Generates harmonized, three-field genotype calls.
- Outputs a consensus genotype combining results from all callers.
- Designed for hg38 reference genome.

#### Inputs
**Runtime Parameters:**
- `String gatk_docker` – Docker image containing GATK
- `String hlahd_docker` – Docker image containing HLA-HD
- `String polysolver_docker` – Docker image containing Polysolver
- `String optitype_docker` – Docker image containing OptiType

**Reference and Input Files:**
- `File original_bam` – Input BAM or CRAM file
- `File original_bam_idx` – Index for the input BAM/CRAM
- `File ref_fasta` – Reference FASTA file
- `File ref_fai` – FASTA index file
- `File ref_dict` – Reference dictionary file
- `File hla_intervals` – Interval list specifying HLA region of chromosome 6

**Support Scripts:**
- `File convert_alleles_python_script` – Script to convert HLA-HD output to three-field calls
- `File count_two_field_alleles_python_script` – Script to count number of two-field calls
- `File hla_groups_file` – File mapping allele names to groups for conversion. This file was obtained from https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom_g.txt.

**Execution Environment Inputs:**
- `String? gcs_project_for_requester_pays` – Optional input that is only required if input files are in a requester-pays bucket
- `File? EMPTY_STRING_HACK` – Placeholder necessary to circumvent the outdated version of cromwell that Terra is tied to. Do not set this input, please leave it as is.

#### Step 1. MakeHLAOnlyBamsAndFastqs
- Extracts reads from HLA regions using GATK.
- Sorts and indexes resulting HLA-only BAM.
- Converts the HLA BAM into paired FASTQs using `samtools fastq`.

#### Step 2. HLAHD
- Runs HLA-HD on paired FASTQ files.
- Post-processes raw results:
  - Standardizes formatting.
  - Converts alleles to three-field resolution.
  - Counts how many alleles are typed to two-field resolution.

#### Step 3. Polysolver (conditional)
- If HLA-HD emits at least one two-field genotype, runs Polysolver. Note: Polysolver only genotypes the HLA-A, HLA-B, and HLA-C genes.
- Converts to standardized format and filters for 3-field resolution.

#### Step 4. Optitype (conditional)
- If HLA-HD emits at least one two-field genotype, runs Optitype.
- Genotypes A/B/C genes using paired-end FASTQs.
- Converts output into a standard format with two alleles per gene.

#### Step 5. Consensus (conditional)
- If both Polysolver and Optitype ran, generates a consensus genotype:
  - For HLA-A, HLA-B, and HLA-C genes, overrides HLA-HD with Polysolver results if the two-field Optitype genotype is consistent with the three-field Polysolver genotype.
  - For other genes, retains HLA-HD results.
  - Produces a final genotype file sorted by gene.
- Note: HLA-HD and Polysolver emit three-field resolution; Optitype emits two-field resolution.

#### Outputs
- `File hlahd_raw_result` – Raw HLA-HD output with reformatted alleles
- `File hlahd_converted_result` – Converted HLA-HD output with group-standardized alleles
- `Int hlahd_two_field_count` – Number of loci with two-field resolution in HLA-HD output
- `Int hlahd_overriden` – Count of loci where HLA-HD was overridden by consensus
- `File consensus` – Final consensus genotype call
- `File? optitype_result` – Optitype genotype result (Optional)
- `File? polysolver_result` – Polysolver genotype result (Optional)

## MakeTable
#### Background

This WDL workflow consolidates HLA consensus calls from multiple samples into a single summary table. The output is a tab-delimited file where each row represents a sample, and each gene is represented by two consecutive columns—one for each allele in the diploid genotype. Homozygous genotypes repeat the same allele in both columns.

Key characteristics:
- Expects precomputed HLA consensus files per sample.
- Sample IDs are provided explicitly in the same order as the consensus files.
- Validates consistency in header format across all samples.
- Produces a TSV-formatted summary table with sample IDs and allele calls.

#### Inputs
- `Array[File] consensus_calls` – List of consensus result files, one per sample
- `Array[String] sample_ids` – List of sample IDs (must be the same length and order as `consensus_calls`)

#### Step 1. Combine
- Parses each consensus file to extract HLA genotype fields.
- Constructs a unified header by inspecting the first sample.
- Appends sample genotype data as rows under the unified header.
- Merges sample IDs into the first column of the final table.

#### Outputs
- `File result` – Final combined summary table (`result.txt`)
