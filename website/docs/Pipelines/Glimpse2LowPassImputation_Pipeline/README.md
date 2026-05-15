---
sidebar_position: 1
slug: /Pipelines/Glimpse2LowpassImputation_Pipeline/README
---

# Imputation Overview

|                                                             Pipeline Version                                                              | Date Updated |        Documentation Author        |                             Questions or Feedback                              |
|:-----------------------------------------------------------------------------------------------------------------------------------------:|:------------:|:----------------------------------:|:------------------------------------------------------------------------------:|
| [Glimpse2LowPassImputation_v0.0.10 (pre-relase)](https://github.com/broadinstitute/warp/releases?q=Glimpse&expanded=true) |  May, 2026   | Terra Scientific Pipeline Services | Please [file an issue in WARP](https://github.com/broadinstitute/warp/issues). |

## Introduction to the GLIMPSE2 Low Pass Imputation pipeline
The GLIMPSE2 Low Pass Imputation pipeline imputes missing genotypes from a list of low-pass CRAM/CRAI files (or a sample manifest pointing to GCS file paths) using a large genomic reference panel. It uses GLIMPSE2 as the imputation tool. Overall, the pipeline splits samples into batches, performs variant calling and imputation on each batch across genomic chunks, and merges the results into a final multi-sample VCF. It outputs the imputed VCF along with key imputation metrics.

## Set-up
 
### Workflow installation and requirements

The [Array Imputation workflow](https://github.com/broadinstitute/warp/blob/develop/pipelines/wdl/arrays/imputation_beagle/ImputationBeagle.wdl) is written in the Workflow Description Language (WDL) and can be deployed using a WDL-compatible execution engine like [Cromwell](https://github.com/broadinstitute/cromwell), a GA4GH compliant, flexible workflow management system that supports multiple computing platforms.
To identify the latest workflow version and release notes, please see the Imputation workflow [changelog](https://github.com/broadinstitute/warp/blob/develop/pipelines/wdl/arrays/imputation_beagle/ImputationBeagle.changelog.md).
The latest release of the workflow, example data, and dependencies are available from the WARP releases page. To discover and search releases, use the WARP command-line tool [Wreleaser](https://github.com/broadinstitute/warp/tree/develop/wreleaser).

Using the Array Imputation pipeline
This pipeline is used by the [All of Us + AnVIL Imputation Service](http://allofus-anvil-imputation.terra.bio). If you choose to use this service, you can impute your samples against the 515,000+ genomes in the All of Us + AnVIL reference panel, which can provide greater accuracy at more sites.

:::tip Try the Imputation pipeline in Terra
You can alternatively run the pipeline with your own panel, using this [WDL](https://github.com/broadinstitute/warp/blob/develop/pipelines/wdl/arrays/imputation_beagle/ImputationBeagle.wdl).
:::
 
### Input descriptions
The table below describes each of the GLIMPSE2 Low Pass Imputation pipeline inputs. The workflow requires a list of .cram, .crai, and sample names (either directly provided or as a manifest). These samples must be from the same species and genotyping chip.

| Input Name                       | Description                                                                                                         |
|----------------------------------|---------------------------------------------------------------------------------------------------------------------|
| `cram_manifest`                  | Optional manifest TSV file containing columns (including header line) of sample_id, cram_path, and cram_index_path referring to cloud-hosted input files to be imputed. This or all three array inputs (crams, cram_indices, sample_ids) must be provided.                               |
| `crams`                          | Optional array of input CRAMs                                                                                       |
| `cram_indices`                   | Optional array of CRAI files corresponding to `crams`                                                               |
| `sample_ids`                     | Optional array of sample ID strings corresponding to CRAM inputs                                                    |
| `contigs`                        | Array of contigs/chromosomes to process                                                                             |
| `reference_panel_prefix`         | Directory/prefix containing `sites.<contig>.vcf.gz`, `sites_table.<contig>.gz`, and `reference_chunks.<contig>.txt` |
| `fasta`                          | Reference FASTA                                                                                                     |
| `fasta_index`                    | FASTA index                                                                                                         |
| `output_basename`                | Basename for intermediate and final outputs                                                                         |
| `ref_dict`                       | Reference dictionary used during ligation/header normalization                                                      |
| `impute_reference_only_variants` | Whether to impute reference-only variants (default: `false`)                                                        |
| `call_indels`                    | Whether to include indels during calling/imputation (default: `false`)                                              |
| `calling_batch_size`             | Batch size for CRAM calling inside each batch subworkflow (default: `100`)                                          |
| `sample_batch_size`              | Batch size at gateway level for splitting very large cohorts (default: `1000`)                                      |
| `gatk_docker`                    | GATK Docker image                                                                                                   |
| `glimpse_docker`                 | GLIMPSE2 Docker image                                                                                               |
| `docker_merge`                   | Docker used for merge/re-annotation step                                                                            |
| `mem_gb_merge`                   | Memory (GB) for post-batch merge/re-annotation (default: `32`)                                                      |

 
## Workflow tasks and tools

The [Array Imputation workflow](https://github.com/broadinstitute/warp/blob/develop/pipelines/wdl/arrays/imputation_beagle/ImputationBeagle.wdl) imports a series of tasks from the ImputationTasks WDL and ImputationBeagleTasks WDL, which are hosted in the Broad [tasks library](https://github.com/broadinstitute/warp/tree/develop/tasks/wdl). The table below describes each workflow task, including the task name, tools, relevant software and non-default parameters. 

| Task / Call                                          | Purpose                                                             | Input Dependencies                                         | Key Function                                                |
|------------------------------------------------------|---------------------------------------------------------------------|------------------------------------------------------------|-------------------------------------------------------------|
| `ConvertCramManifestToInputArrays`                   | Convert cram manifest input into CRAMs/CRAIs/sample IDs arrays      | `cram_manifest` | Facilitates submission of very large sets of inputs via manifest file           |
| `SplitIntoSampleBatches`                             | Split CRAMs/CRAIs/sample IDs into sample-level batches              | `crams`, `cram_indices`, `sample_ids` from inputs or derived from `cram_manifest`, `sample_batch_size` | Enables large-cohort scaling at gateway level               |
| `RunBatch` (`Glimpse2LowPassImputationBatch`)        | Run full low-pass imputation pipeline on each sample batch          | Batch-specific CRAMs/indices/sample IDs + reference inputs | Produces per-batch, per-contig ligated imputed VCFs         |
| `ExtractAnnotations`                                 | Extract AF/INFO annotations from each batch contig VCF              | Batch ligated VCFs and indexes                             | Captures annotations needed for post-merge recomputation    |
| `MergeContigVcfs` (`MergeSampleChunksVcfsWithPaste`) | Merge sample columns across batch VCFs for one contig               | Array of batch VCFs for contig                             | Creates full-cohort contig VCF with aligned site lists      |
| `RecomputeAndAnnotate`                               | Recompute AF/INFO across merged cohort and write updated contig VCF | Merged contig VCF + extracted annotations                  | Restores cohort-correct annotations after paste-based merge |
| `SelectContigVariants`                               | Create variants-only contig VCF                                     | Re-annotated contig VCF                                    | Removes homozygous-reference-only records                   |
| `CreateContigHomRefVcf`                              | Create hom-ref-sites-only contig VCF                                | Re-annotated contig VCF                                    | Keeps homozygous-reference-only sites                       |
| `CombineBatchCoverageMetrics`                        | Combine optional coverage metric files across batches               | `RunBatch.coverage_metrics`                                | Produces aggregated coverage table when metrics exist       |
| `GatherVcfsNoIndex`                                  | Gather contig variant VCFs into genome-wide variant VCF             | Variant-only contig VCFs                                   | Produces final genome-wide variant VCF                      |
| `CreateVcfIndexAndMd5`                               | Index and checksum final variant VCF                                | Gathered variant VCF                                       | Creates `.tbi` and md5                                      |
| `GatherVcfsNoIndexHomRefOnly`                        | Gather contig hom-ref-sites-only VCFs                               | Hom-ref contig VCFs                                        | Produces final genome-wide hom-ref-sites-only VCF           |
| `CreateVcfIndexAndMd5HomRefOnly`                     | Index and checksum final hom-ref-sites-only VCF                     | Gathered hom-ref-sites-only VCF                            | Creates `.tbi` and md5                                      |
| `CollectQCMetrics`                                   | Compute sample QC metrics from final imputed variant VCF            | Final imputed variant VCF                                  | Generates sample-level QC report                            |


## Workflow outputs

The table below summarizes the workflow outputs. If running the workflow on Cromwell, these outputs are found in the task execution directory.
 
| Output name                          | Description                                                                                                                                                                                                                  | Type  |
|--------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------|
| imputed_multi_sample_vcf             | Final imputed multi-sample WGS VCF file, only contains variant records                                                                                                                                                       |
| imputed_multisample_vcf_index        | Index file for VCF from the CreateIndexForGatheredVcf task.                                                                                                                                                                  | Index |
| imputed_hom_ref_sites_only_vcf       | Sites only VCF that contains sites where all samples are homozygous reference                                                                                                                                                |
| imputed_hom_ref_sites_only_vcf_index | Index file for hom ref sites only VCF                                                                                                                                                                                        |
| chunks_info                          | TSV from StoreMetricsInfo task; contains the chunk intervals, variant counts per chunk (filtered input and overlap with reference panel), and whether each chunk was successfully imputed.                                   | TSV   |
| contigs_info                         | TSV from StoreMetricsInfo task; contains contig-level (chromosome-level) aggregated metrics including total variant counts in the raw input, filtered input, and overlap with reference panel for each processed chromosome. | TSV   |


 
## Important notes
 
- Runtime parameters are optimized for Broad's Google Cloud Platform implementation.

## Citing the Imputation Pipeline

If you use the GLIMPSE2 Low Pass Imputation Pipeline in your research, please consider citing our preprint:

Degatano, K., Awdeh, A., Cox III, R.S., Dingman, W., Grant, G., Khajouei, F., Kiernan, E., Konwar, K., Mathews, K.L., Palis, K., et al. Warp Analysis Research Pipelines: Cloud-optimized workflows for biological data processing and reproducible analysis. Bioinformatics 2025; btaf494. https://doi.org/10.1093/bioinformatics/btaf494
 
## Contact us

Help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.
 
## Licensing
 
Copyright Broad Institute, 2020 | BSD-3
 
The workflow script is released under the **WDL open source code license (BSD-3)** (full license text at https://github.com/broadinstitute/warp/blob/master/LICENSE). However, please note that the programs it calls may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running this script.
 
- [GATK](https://github.com/broadinstitute/gatk/blob/master/LICENSE.TXT)
- [Picard](https://github.com/broadinstitute/picard/blob/master/LICENSE.txt)
- [GLIMPSE](https://github.com/odelaneau/GLIMPSE/blob/master/LICENSE)
- [bcftools](https://github.com/samtools/bcftools/blob/develop/LICENSE)
- [vcftools](http://vcftools.sourceforge.net/license.html)
 
 


