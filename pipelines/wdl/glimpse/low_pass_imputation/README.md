## GLIMPSE2 Low-Pass Imputation Summary

The `Glimpse2LowPassImputation` workflow is a WDL-based pipeline for low-pass whole genome imputation using [GLIMPSE2](https://odelaneau.github.io/GLIMPSE/).  
This top-level workflow is now a gateway that scales to large cohorts by splitting samples into batches, running a per-batch imputation subworkflow, then merging batch outputs back into cohort-level results.

The workflow processes each requested contig independently, imputes each sample batch against reference-defined chunks, ligates chunk outputs per batch/contig, merges sample columns across batches, recomputes AF/INFO annotations, and gathers contig outputs into final genome-wide files.

### Pipeline Features

| Pipeline features       | Description                                                                           | Source                                                                   |
|-------------------------|---------------------------------------------------------------------------------------|--------------------------------------------------------------------------|
| Assay type              | Low-pass whole genome imputation using GLIMPSE2                                       | [GLIMPSE2](https://odelaneau.github.io/GLIMPSE/)                         |
| Overall workflow        | CRAM calling, shard-based phasing/imputation, ligation, batch merge, and QC           | Defined in `Glimpse2LowPassImputation.wdl` + imported subworkflows/tasks |
| Workflow language       | WDL 1.0                                                                               | [openWDL](https://github.com/openwdl/wdl)                                |
| Sub-workflows           | Gateway workflow + `Glimpse2LowPassImputationBatch`                                   | Imported from `Glimpse2LowPassImputationBatch.wdl`                       |
| Genomic processing      | Contig-by-contig and reference-chunk-based processing                                 | Workflow scatter logic                                                   |
| Cohort scalability      | Inputs optionally specified via `cram_manifest`, sample batching (`sample_batch_size`) then per-contig batch merge                     | Gateway orchestration in `Glimpse2LowPassImputation.wdl`                 |
| Algorithms              | `bcftools` mpileup/call/norm/merge + GLIMPSE2 phase/ligate + post-merge re-annotation | Task commands in batch and task WDLs                                     |
| Quality control         | Sample QC metrics and optional coverage-metrics aggregation                           | `CollectQCMetrics`, `CombineCoverageMetrics`                             |
| Data input file format  | CRAM/CRAI arrays with sample IDs                                                      | Workflow input block                                                     |
| Data output file format | Imputed VCFs, indexes, md5s, and QC/coverage metric tables                            | Workflow outputs                                                         |
| Containers              | GATK, GLIMPSE2, bcftools/samtools suite, Hail, Python, Ubuntu                         | Runtime blocks                                                           |
| Resource optimization   | Parallelization by sample batch, contig, and reference shard                          | Workflow architecture                                                    |

### Inputs

This gateway workflow expects CRAM-based inputs and a GLIMPSE2-compatible reference panel layout.

| Input                            | Description                                                                                                         |
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

### Workflow Tasks

The top-level workflow orchestrates batching, per-batch imputation, and cohort-level merging/re-annotation.

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

### Outputs

Upon successful completion, the workflow emits final genome-wide imputed outputs, corresponding index and checksum files, and QC metrics. Coverage metrics are optional.

| Output                                 | Description                                                     |
|----------------------------------------|-----------------------------------------------------------------|
| `imputed_vcf`                          | Final imputed multi-sample variant VCF                          |
| `imputed_vcf_index`                    | Index file for final imputed VCF                                |
| `imputed_vcf_md5sum`                   | MD5 checksum for final imputed VCF                              |
| `imputed_hom_ref_sites_only_vcf`       | Final sites-only VCF containing homozygous-reference-only sites |
| `imputed_hom_ref_sites_only_vcf_index` | Index file for hom-ref-sites-only VCF                           |
| `imputed_hom_ref_sites_only_vcf_md5`   | MD5 checksum for hom-ref-sites-only VCF                         |
| `qc_metrics`                           | Sample-level QC metrics table                                   |
| `coverage_metrics`                     | Optional combined coverage metrics table                        |


## Glimpse2LowPassImputationBatch summary

The `Glimpse2LowPassImputationBatch` workflow is the per-batch subworkflow used by the top-level `Glimpse2LowPassImputation` gateway workflow.  
It is designed for cohorts up to roughly 1000 samples per batch, then returns contig-level ligated imputed VCFs that can be merged across batches upstream.

### Batch Workflow Role

- runs low-pass variant calling from CRAMs at reference-panel sites
- phases/imputes each contig in reference-defined chunks using GLIMPSE2
- ligates chunk-level outputs back into one imputed VCF per contig
- emits optional coverage metrics aggregated across chunks/contigs

### Batch Inputs

| Input                            | Description                                                                                               |
|----------------------------------|-----------------------------------------------------------------------------------------------------------|
| `contigs`                        | Contigs/chromosomes to process                                                                            |
| `reference_panel_prefix`         | Prefix containing `sites.<contig>.vcf.gz`, `sites_table.<contig>.gz`, and `reference_chunks.<contig>.txt` |
| `crams`                          | CRAM files for this batch                                                                                 |
| `cram_indices`                   | CRAI files for the batch CRAMs                                                                            |
| `sample_ids`                     | Sample IDs aligned to `crams`                                                                             |
| `fasta` / `fasta_index`          | Reference FASTA and index                                                                                 |
| `output_basename`                | Basename for intermediate and emitted files                                                               |
| `ref_dict`                       | Reference dictionary used during ligation/reheader                                                        |
| `impute_reference_only_variants` | Pass-through option for GLIMPSE2 phase                                                                    |
| `call_indels`                    | Whether to include indels during calling/imputation                                                       |
| `calling_batch_size`             | Internal batch size for CRAM calling fan-out within this subworkflow                                      |
| `gatk_docker` / `glimpse_docker` | Container images for GATK and GLIMPSE2 tools                                                              |

### Batch Internal Processing

| Step                                   | Purpose                                                                               |
|----------------------------------------|---------------------------------------------------------------------------------------|
| `SplitIntoBatches` (conditional)       | Splits CRAMs/CRAIs/sample IDs into internal calling batches                           |
| `BcftoolsMpileup`                      | Computes pileups at panel sites per internal batch                                    |
| `BcftoolsCall`                         | Calls candidate variants from mpileup output                                          |
| `BcftoolsNorm`                         | Normalizes and indexes called variants                                                |
| `BcftoolsMerge` (conditional)          | Merges per-internal-batch VCFs if multiple were produced                              |
| `ComputeShardsAndMemoryPerShard`       | Reads reference chunks and computes per-shard memory estimates                        |
| `GlimpsePhase`                         | Runs `GLIMPSE2_phase` for each reference shard                                        |
| `GlimpseLigate`                        | Ligates shard outputs to one contig-level imputed VCF and updates sequence dictionary |
| `CombineCoverageMetrics` (conditional) | Combines optional shard/contig coverage metric files                                  |

### Batch Outputs

| Output                               | Description                                           |
|--------------------------------------|-------------------------------------------------------|
| `imputed_contig_ligated_vcfs`        | Per-contig ligated imputed VCFs for this sample batch |
| `imputed_contig_ligated_vcf_indices` | Index files for each contig ligated VCF               |
| `coverage_metrics`                   | Optional batch-level combined coverage metrics table  |

### Glimpse2MergeBatches AF and INFO score recalculation
**By: Christopher Kachulis**

Glimpse outputs three values in the info field:
1) AF.  This is easily calculated for the full cohort based on the individual batch values as a weighted mean of the batch allele frequencies, $$AF_{cohort}=\frac{\sum AF_i N_i}{\sum N_i}$$

2) RAF.  This is the reference panel allele frequency.  Assuming the same reference panel was used for all batches, these are all the same, so just take the first value

3) INFO.  This is the "IMPUTE style info score".  Based on https://static-content.springer.com/esm/art%3A10.1038%2Fnrg2796/MediaObjects/41576_2010_BFnrg2796_MOESM3_ESM.pdf, this is calculated as $$1-\frac{\sum f_j - e^2_j}{2N\times AF(1-AF)}$$, or 1 if AF=0,1, with $j$ running across the N samples, $f_j = p_{j1} + 4 p_{j2}$, $e_j = p_{j1} + 2 p_{j2}$, where $p_{j1}$ is the imputed posterior of a het for sample $j$, and $p_{j2}$ is the imputed posterior of a hom var for sample $j$.  Note that the terms in the denominator are all easily calculated for a cohort based on their values for the constituent batches;  AF as described above, and N just as the sum over the batches.  The numerator we can define for batch $i$ as $C_i = \sum (f_{ij} -e_{ij}^2)$ and note that, for a whole cohort, we simply have $C_{cohort} = \sum C_i$.  We then note that, for a batch i, we have $$I_i =1 - \frac{C_i}{2N_i\times AF_i(1-AF_i)}$$, so we can solve for $C_i$ as $$C_i = (1-I_i)*2N_i\times AF_i(1-AF_i)$$.  We can then calculate the cohort INFO score as $$I_{cohort}=1-\frac{\sum C_i}{2N_{cohort} AF_{cohort}(1-AF_{cohort})}$$ which becomes, $$I_{cohort}=1-\frac{\sum (1-I_i)*2N_i\times AF_i(1-AF_i)}{2\sum N_i \times \frac{\sum AF_i N_i}{\sum N_i} (1-\frac{\sum AF_i N_i}{\sum N_i})}$$


## Glimpse2LowPassImputationQuotaConsumed summary

The `QuotaConsumed` workflow computes submitted sample count for service quota accounting.  
Quota is derived from number of CRAM entries found in the provided CRAM manifest.

## Glimpse2LowPassImputationQC summary

The `InputQC` workflow validates CRAM-based inputs supplied by the `cram_manifest` for GLIMPSE2 low-pass imputation. Checks include:

- required manifest columns are present: `sample_id`, `cram_path`, `cram_index_path`
- counts of CRAMs, CRAIs, and sample IDs match
- sample IDs are unique
- CRAM paths are unique
- CRAM file names end with `.cram` and indices end with `.crai`
- input paths use `gs://` format and are accessible
- CRAM file sizes do not exceed the configured maximum (default: 10 GB)
- optional requester-pays validation via `billing_project_for_rp`
