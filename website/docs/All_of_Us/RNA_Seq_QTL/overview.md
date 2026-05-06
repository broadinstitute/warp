---
title: All of Us RNA-seq eQTL and sQTL Analysis Pipeline
sidebar_position: 1
className: aou-doc-page
---

<div className="aou-folder-text">

This section documents the All of Us RNA-seq QTL workflows for expression/splicing phenotype generation, QTL analysis preparation, and SuSiE fine-mapping aggregation.

## Quick Summary

* **Purpose:** Run reproducible eQTL and sQTL analyses from genotype and RNA-derived inputs.
* **Primary Outputs:** RNA-derived phenotype artifacts, TensorQTL-ready inputs, SuSiE fine-mapping results, and aggregated annotations.

## Input Requirements

To run the workflows described below, you generally need:

* A joint-called **VCF** containing the relevant samples
* **Research IDs** partitioned by ancestry or subpopulation
* RNA expression quantifications (for eQTL)
* BAM/CRAM files for splice junction extraction (for sQTL)
* Sample-level metadata tables

## Ordered Analysis Flow

The table below reflects the original end-to-end run order, including steps that do not yet have dedicated docs pages in this folder.

| Order | Stage | Workflow / Component | Documentation | WDL | Run next |
| :--: | --- | --- | --- | --- | --- |
| 0 | Cohort setup | Ancestry grouping and sample lists | No dedicated WDL page (notebook/table prep) | N/A | Use ancestry/sample partitions as inputs to genotype and phenotype prep. |
| 1 | Genotype prep | Prepare genotypes (pruning + PLINK + PCs) | No dedicated page yet | [PrepareGenotypes.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/PrepareGenotypes.wdl) | Run dosage generation per ancestry/population. |
| 2 | Genotype prep | Calculate genotype dosage | No dedicated page yet | [calculateGenotypeDosage.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/prepare_QTL/calculateGenotypeDosage.wdl) | Feed dosages into TensorQTL and SuSiE inputs later. |
| 3 | RNA processing | RNA-seq AoU processing (alignment/quant/QC) | [RNA-seq AoU Processing](./rnaseq_aou) | [rnaseq_aou.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/GTEx/rnaseq_aou.wdl) | Branch into eQTL phenotype prep and/or sQTL junction extraction. |
| 4A | eQTL phenotypes | Prepare eQTL phenotype BED + phenotype PCs | No dedicated page yet | [prepare_eQTL.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/prepare_QTL/prepare_eQTL.wdl) | Merge covariates for eQTL TensorQTL run. |
| 4B | sQTL phenotypes | Extract junctions from BAM | [Leafcutter BAM to Junctions](./leafcutter_bam_to_junc) | [leafcutter_bam_to_junc.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/leafcutter_bam_to_junc.wdl) | Cluster junctions to build sQTL phenotype matrices. |
| 5B | sQTL phenotypes | Cluster junctions + generate leafcutter outputs | [Leafcutter Clustering](./leafcutter_cluster) | [leafcutter_cluster.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/leafcutter_cluster.wdl) | If needed, run separate phenotype-group generation before covariate merge. |
| 5C | sQTL phenotypes | Prepare sQTL phenotype BED + PCs | No dedicated page yet | [prepare_sQTL.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/prepare_QTL/prepare_sQTL.wdl) | Use splicing BED/PC outputs for sQTL covariate merge and TensorQTL. |
| 6B | sQTL metadata | Calculate phenotype groups | [Calculate Phenotype Groups](./calculate_phenotype_groups) | [CalculatePhenotypeGroups.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/CalculatePhenotypeGroups.wdl) | Merge covariates for sQTL TensorQTL run. |
| 7 | Covariates | Merge covariates (genotype PCs + phenotype PCs ± groups) | No dedicated page yet | [MergeCovariates.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/prepare_QTL/MergeCovariates.wdl) | Run TensorQTL cis permutations for eQTL/sQTL. |
| 8 | Association | TensorQTL cis permutations | No dedicated page yet | [tensorqtl_cis_permutations.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/tensorQTL_cis_permutations/tensorqtl_cis_permutations.wdl) | Recalculate FDR and prepare significant loci for fine-mapping. |
| 9 | Fine-mapping prep | FDR recalculation + SuSiE input preparation | No dedicated WDL page in this folder | N/A | Run SuSiE per phenotype window. |
| 10 | Fine-mapping | SuSiE fine-mapping | [SuSiE Fine-Mapping Workflow](./susieR_workflow) | [susieR_workflow.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/susieR_workflow.wdl) | Aggregate SuSiE outputs across phenotypes. |
| 11 | Aggregation | Aggregate SuSiE outputs and annotate | [Aggregate SuSiE Workflow](./aggregate_susie_workflow) | [AggregateSusieWorkflow.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/AggregateSusieWorkflow.wdl) | Optional: compute AF summaries for downstream interpretation/reporting. |
| 12 (optional) | Variant summaries | Calculate allele frequencies | No dedicated page yet | [calculateAF.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/prepare_QTL/calculateAF.wdl) | Use AF metrics in interpretation/QC/reporting. |

## Practical Run Notes

* **eQTL path:** 0 → 1 → 2 → 3 → 4A → 7 → 8 → 9 → 10 → 11 (+12 optional).
* **sQTL path:** 0 → 1 → 2 → 3 → 4B → 5B → 5C → 6B (if required) → 7 → 8 → 9 → 10 → 11 (+12 optional).
* **Key dependency:** `susieR_workflow` expects TensorQTL-derived significant loci plus dosage inputs from earlier genotype steps.
* **Phenotype groups:** required for many sQTL TensorQTL configurations; for eQTL they are typically not required.

## Additional Processing Notes

* **FDR and SuSiE prep:** after TensorQTL, recalculate FDR, filter significant loci (commonly 0.05), and format SuSiE-ready inputs.
* **SuSiE runtime guidance:** preemptible VMs can reduce cost; pinned Docker SHAs improve reproducibility.
* **Aggregation inputs:** for `AggregateSusieWorkflow`, use fine-mapped `SusieParquet` outputs; do not use full/all-tested parquet outputs.

## Acknowledgements

The original versions of these workflows were either created by the GTEx Consortium (see their [GTEx GitHub repository](https://github.com/broadinstitute/gtex-pipeline/tree/master?tab=readme-ov-file)) or by the **Dr. Stephen Montgomery Lab** at Stanford University. Most of the eQTL scripts originated from the [AoU-Multiomics-Analysis repository](https://github.com/AoU-Multiomics-Analysis).

This pipeline builds upon extensive work by the **Stephen Montgomery Lab** at Stanford University. Special thanks to:

* **Evin Padhi**
* **Jon Nguyen**

for developing foundational versions of many scripts and workflows used in this analysis.

Additional integration, optimization, and workflow migration were performed by the All of Us Multiomics and the Broad Pipeline Development teams as part of the WARP workflow suite.

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>

