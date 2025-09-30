

## Imputation Beagle Summary

The ImputationBeagle workflow (v2.0.0) is a WDL-based pipeline for genotype imputation using [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html). This workflow processes multi-sample VCF files, performing phasing and imputation by chromosome chunks to increase genetic variant resolution.

The workflow divides each chromosome into chunks, phases and imputes variants using the reference panel, and then combines results back into a complete multi-sample VCF.

### Pipeline Features

This workflow leverages Beagle's advanced algorithms for phasing and imputation while incorporating robust quality control measures and scalable processing architecture.

| Pipeline features       | Description                                                                                         | Source                                                                        |
|-------------------------|-----------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------|
| Assay type              | Genotype imputation using phased reference panels                                                   | [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html)          |
| Overall workflow        | Sample counting, chunking, phasing, imputation, and post-processing                                 | Code available in the provided WDL script                                     |
| Workflow language       | WDL 1.0                                                                                             | [openWDL](https://github.com/openwdl/wdl)                                     |
| Sub-workflows           | Imports tasks from ImputationTasks.wdl and ImputationBeagleTasks.wdl                                | Referenced in the workflow                                                    |
| Genomic processing      | Chromosome-based chunking with overlap to reduce edge effects                                       | Defined in workflow processing logic                                          |
| Reference resources     | Reference panels (bref3 format), genetic maps, and reference dictionary                             | Referenced as input parameters                                                |
| Algorithms              | Beagle phasing and imputation                                                                       | [Browning et al.](https://faculty.washington.edu/browning/beagle/beagle.html) |
| Quality control         | Variant counting and chunk validation before costly operations                                      | Implemented in CheckChunks task                                               |
| Scalability             | Dynamic CPU/memory allocation based on sample count                                                 | Resource scaling logic in workflow                                            |
| Data input file format  | Multi-sample Array VCF                                                                              | Standard genotype format                                                      |
| Data output file format | Imputed VCF with index and QC information                                                           | Standard output formats                                                       |
| Post-processing         | Updating headers                                                                                    | GATK-based processing                                                         |
| Containers              | GATK Docker (4.6.0.0) and Ubuntu Docker (20.04)                                                     | Specified in the workflow                                                     |
| Resource optimization   | Parallelization by chromosome and chunks both positionally and by sample with preemptible execution | Architecture of workflow                                                      |

### Inputs

The Beagle workflow requires several essential inputs to perform genotype imputation effectively. These inputs include the source VCF data, reference materials for alignment and phasing, and configuration parameters that control the chunking strategy and processing behavior.

| Input                         | Description                                                              |
|-------------------------------|--------------------------------------------------------------------------|
| `multi_sample_vcf`            | Multi-sample VCF file containing genotype data                           |
| `ref_dict`                    | Reference dictionary for contig information and header updating          |
| `contigs`                     | Array of allowed contigs/chromosomes to process                          |
| `reference_panel_path_prefix` | Path and prefix to reference panel files in bucket                       |
| `genetic_maps_path`           | Path to genetic maps for all contigs                                     |
| `output_basename`             | Basename for intermediate and output files                               |
| `chunkLength`                 | Length of chunks for processing (default: 25,000,000)                    |
| `chunkOverlaps`               | Overlap padding to reduce edge effects (default: 2,000,000)              |
| `sample_chunk_size`           | Number of samples to chunk by when processing (default: 1,000)           |
| `bed_suffix`                  | File extension for reference panel bed files (default: .bed)             |
| `bref3_suffix`                | File extension for reference panel bref3 files (default: .bref3)         |
| `beagle_cpu`                  | Number of cpus to use for Beagle Phase and Impute tasks (default: 8)     |
| `beagle_phase_memory_in_gb`   | Memory, in GB, to use for Beagle Phase task (default: 40)                |
| `beagle_impute_memory_in_gb`  | Memory, in GB, to use for Beagle Impute task (default: 45)               |
| `gatk_docker`                 | GATK Docker image (default: us.gcr.io/broad-gatk/gatk:4.6.0.0)           |
| `ubuntu_docker`               | Ubuntu Docker image (default: us.gcr.io/broad-dsde-methods/ubuntu:20.04) |
| `error_count_override`        | Optional override for error count threshold                              |

### Workflow Tasks

The Beagle workflow consists of multiple interconnected tasks that work together to process VCF data through preparation, quality control, imputation, and post-processing stages.

| Task                                   | Purpose                                                                            | Input Dependencies                              | Key Function                                                                 |
|----------------------------------------|------------------------------------------------------------------------------------|-------------------------------------------------|------------------------------------------------------------------------------|
| `CountSamples`                         | Count number of samples in input VCF                                               | multi_sample_vcf                                | Determines resource allocation for downstream tasks                          |
| `CreateVcfIndex`                       | Index the input VCF file                                                           | multi_sample_vcf                                | Creates index for efficient VCF access                                       |
| `CalculateContigsToProcess`            | Determine which contigs will be processed by the workflow                          | multi_sample_vcf, contigs                               | Extracts contigs from input VCF and filters by allowed contigs                     |
| `CalculateChromosomeLength`            | Calculate length of each chromosome                                                | ref_dict, contig                                | Determines number of chunks needed per chromosome                            |
| `GenerateChunk`                        | Create chunked VCF files with overlaps                                             | indexed VCF, coordinates                        | Splits chromosomes into processable chunks                                   |
| `CountVariantsInChunks`                | Count variants in chunks vs reference panel                                        | chunk VCF, reference panel                      | Quality control for chunk validity                                           |
| `CheckChunks`                          | Validate chunk quality                                                             | variant counts                                  | Ensures chunks meet quality thresholds                                       |
| `StoreChunksInfo`                      | Store chunk metadata and QC metrics                                                | chunk information                               | Tracks processing statistics                                                 |
| `ErrorWithMessageIfErrorCountNotZero`  | Fail workflow if CheckChunks fail                                                  | error counts                                    | Quality control gate before expensive operations                             |
| `SelectSamplesWithCut`                 | Splits vcf into chunks of samples using `cut`                                      | chunk VCF                                       | helps with scaling the wdl                                                   |
| `Phase`                                | Phase genotypes using Beagle                                                       | chunk VCF, reference panel, genetic map         | Determines haplotype phase for imputation                                    |
| `Impute`                               | Impute missing genotypes using Beagle                                              | phased VCF, reference panel                     | Main imputation step                                                         |
| `LocalizeAndSubsetVcfToRegion`         | Remove overlap region from imputed VCF                                             | imputed VCF                                     | Makes sure variants are not counted multiple times due to padding            |
| `QuerySampleChunkedVcfForReannotation` | Retrieve DS/AP1/AP2 annotations from sample chunked VCFs for processing AF and DR2 | imputed VCF                                     | Part of the merging process of sample chunked VCFs                           |
| `RemoveAPAnnotations`                  | Remove AP1 and AP2 annotations from imputed VCF                                    | imputed VCF                                     | Cuts down on size of vcf for following intermediate steps                    |
| `RecalculateDR2AndAFChunked`           | Recalculates DR2 and AF for sample merged VCF                                      | QueryMergedVcfForReannotation output query file | Summarizes statistics necessary to merge back sample chunked VCFs            |
| `MergeSampleChunksVcfsWithPaste`       | Merge sample chunked VCFs with `Paste`                                             | sample chunked VCFs                             | Merges sample chunked VCFs back together                                     |
| `IndexMergedSampleChunksVcfs`          | Create an index for merged sample merged VCF                                       | merged sample chunked VCFs                      | Indexes are useful for downstream processing                                 |
| `AggregateChunkedDR2AndAF`             | Calculate site level values for AF and DR2 using split chunk numbers               | sample chunk DS, AP1, AP2 summarized values     | Allows us to calculate site level annotations across the sample chunked vcfs |
| `ReannotateDR2AndAF`                   | apply calculated AF and DR2 annotations to merged VCF                              | Aggregated DR2 and AF annotation file           | We now have a fully merged, correctly annotated vcf                          |
| `UpdateHeader`                         | Update VCF header with reference info                                              | imputed VCF, ref_dict                           | Ensures proper VCF formatting                                                |
| `GatherVcfsNoIndex`                    | Combine all chromosome chunks                                                      | all processed chunks                            | Creates final output VCF                                                     |
| `IndexMergedSampleChunksVcfs`          | Create an index for final output                                                   | final output VCf                                | Indexes are useful for downstream processing                                 |


### Outputs

Upon successful completion, the workflow produces a fully imputed multi-sample VCF file along with supporting files for downstream analysis. The primary outputs include the final imputed VCF with its index file and some quality control information about the processing chunks.

| Output                           | Description                                          |
|----------------------------------|------------------------------------------------------|
| `imputed_multi_sample_vcf`       | Final imputed multi-sample WGS VCF file              |
| `imputed_multi_sample_vcf_index` | Index file for the imputed VCF                       |
| `chunks_info`                    | Information about processed chunks and QC statistics |


## ArrayImputationQuotaConsumed summary

The ArrayImputationQuotaConsumed pipeline is used by the All of Us/AnVIL Imputation Service and calculates the number of samples in the input multi-sample VCF, which is the metric used by the service for ImputationBeagle pipeline quota.


## ArrayImputationQC summary

The ArrayImputationQuotaConsumed pipeline is used by the All of Us/AnVIL Imputation Service and runs various qc checks on the input multi-sample VCF, checks include:
- vcf version 4.x
- hg38 chromosome names
- variants on each of the hg38 canonical chromosomes.
- is not a WGS array vcf
