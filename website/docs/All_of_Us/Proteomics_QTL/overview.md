---
title: All of Us Proteomics pQTL Analysis Pipeline
sidebar_position: 1
className: aou-doc-page
---

<div className="aou-folder-text">

:::caution WORK IN PROGRESS
The All of Us proteomics pQTL workflows and documentation are still in progress. File locations, workflow interfaces, container versions, and analysis parameters may change as the process is finalized.
:::

This section documents the All of Us (AoU) proteomics pQTL analysis process. The workflow uses AoU genetic ancestry assignments to create ancestry-specific and combined analysis inputs. Several genotype-preparation and TensorQTL steps are shared with the [RNA-seq QTL workflows](../RNA_Seq_QTL/overview).

The foundational analyses were developed by the [Dr. Stephen Montgomery Lab](https://med.stanford.edu/montgomerylab.html) at Stanford University and are based in part on work from the [AoU Multiomics Analysis repository](https://github.com/AoU-Multiomics-Analysis). [See acknowledgements](#acknowledgements).

## Quick Summary

* **Purpose:** Prepare proteomic phenotypes and covariates and run ancestry-specific or combined cis-pQTL analyses.
* **Primary inputs:** AoU genotype data, genetic ancestry assignments, a reference-plate median-normalized and replicate-removed proteomics matrix, and genetic principal components.
* **Primary outputs:** Proteomics phenotype BED files, proteomic phenotype PCs, merged QTL covariates, and TensorQTL cis-association results.
* **Current status:** In progress; this page describes the current analysis process rather than a finalized production pipeline.

## Input Requirements

The current process generally requires:

* An AoU proteomics release with normalized NPX measurements
* A reference-plate median-normalized matrix with technical replicates removed, such as `rep_removed.tsv`, containing all analysis samples
* Updated normalized values in the `NPX` columns and normalization adjustment factors in the `adj_factor` column
* AoU genetic ancestry assignments for the analysis populations
* An AoU genotype MatrixTable or exported genotype VCF
* A participant sample list shared with the proteomics data
* Genetic principal components and PLINK2 genotype files for each analysis population
* GENCODE v48 GRCh38 gene annotations for protein-to-gene coordinate mapping
* Versioned cloud-storage locations for intermediate and final outputs

The current population labels are `AFR`, `AMR`, `EAS`, `EUR`, `MID`, `SAS`, and `COMB` for the combined analysis. The exact population definitions and sample counts should be recorded with each analysis because they depend on the AoU genetic ancestry assignments and input release.

## Ordered Analysis Flow

The table below describes the current end-to-end pQTL preparation, association, and fine-mapping flow. Steps marked as shared use workflows that are also used by the AoU RNA-seq QTL process.

| Order | Stage | Workflow / Component | WDL or source | Run next |
| :--: | --- | --- | --- | --- |
| 0 | Cohort setup | Apply the reference-plate median approach, remove technical replicates, and create the all-sample `rep_removed.tsv` matrix. Use the updated normalized values in the `NPX` columns and retain normalization adjustment factors in `adj_factor`; then select the analysis cohort and partition samples by AoU genetic ancestry. | No cohort-setup script is currently available in WARP. | Provide the replicate-removed matrix and consistent ancestry-specific sample lists to genotype and phenotype preparation. |
| 1 | Genotype prep | Filter the AoU MatrixTable, export a VCF, create PLINK2 files, and calculate genetic PCs | [PrepareGenotypes.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/PrepareGenotypes.wdl) (shared with RNA-seq QTL) | Prepare dosage and genotype inputs for association testing. |
| 2 | Genotype prep | Calculate genotype dosage for the analysis cohort | [calculateGenotypeDosage.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/prepare_QTL/calculateGenotypeDosage.wdl) (shared with RNA-seq QTL) | Use the dosage files with TensorQTL and downstream analyses. |
| 3 | Proteomics phenotypes | Filter complete samples, identify proteomic outliers, map UniProt identifiers to gene coordinates, and build ancestry-specific phenotype matrices | [Proteomics_prepare_covariates.ipynb](https://github.com/broadinstitute/warp/blob/develop/all_of_us/proteomics/Proteomics_prepare_covariates.ipynb) | Continue with phenotype normalization and covariate preparation. |
| 4 | Proteomics phenotypes | Rank-normalize protein measurements within each ancestry group and calculate proteomic phenotype PCs | [Proteomics_prepare_covariates.ipynb](https://github.com/broadinstitute/warp/blob/develop/all_of_us/proteomics/Proteomics_prepare_covariates.ipynb) | Merge phenotype PCs with genetic PCs. |
| 5 | Covariates | Generate proteomic phenotype PCs, combine them with genetic PCs, and add the resulting covariate paths to the mapping inputs | [Proteomics_prepare_covariates.ipynb](https://github.com/broadinstitute/warp/blob/develop/all_of_us/proteomics/Proteomics_prepare_covariates.ipynb); mapping-input table updates are manual | Supply the merged covariates to TensorQTL. |
| 6 | Association | Run cis-TensorQTL permutations using PLINK2 genotype inputs, proteomic phenotype BED files, and merged covariates | [tensorqtl_cis_permutations.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/tensorQTL_cis_permutations/tensorqtl_cis_permutations.wdl) (shared with RNA-seq QTL) | Review TensorQTL results and apply the planned significance and quality-control filters. |
| 7 | Association QC | Compare cis-pQTL results across releases and ancestries and perform combined pQTL quality control | Not yet represented by a dedicated WDL in WARP. | Select results for downstream fine-mapping and reporting. |
| 8 | Fine-mapping preparation | Create and update the mapping-input table, prepare phenotype/genotype/variant metadata, and configure the fine-mapping inputs | Mapping-input preparation is manual and fine-mapping preparation is in progress; no dedicated WDL is currently available in WARP | Run SuSiE by ancestry and for the combined cohort. |
| 9 | Fine-mapping | Run SuSiE fine-mapping for each ancestry and the combined analysis | [susieR_workflow.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/susieR_workflow.wdl) (shared with RNA-seq QTL) | Calculate allele frequencies and aggregate fine-mapping outputs. |
| 10 | Fine-mapping annotation | Calculate allele frequencies for fine-mapping and annotation | [calculateAF.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/prepare_QTL/calculateAF.wdl) (shared with RNA-seq QTL) | Provide allele-frequency outputs to SuSiE aggregation. |
| 11 | Aggregation | Aggregate annotated SuSiE results across phenotypes and populations | [AggregateSusieWorkflow.wdl](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/AggregateSusieWorkflow.wdl) (shared with RNA-seq QTL) | Run combined pQTL and SuSiE QC and prepare public analysis outputs. |

## Proteomics Phenotype and Covariate Preparation

The current proteomics-specific preparation is implemented in [`Proteomics_prepare_covariates.ipynb`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/proteomics/Proteomics_prepare_covariates.ipynb). The notebook is an R notebook despite its `.ipynb` format and assumes that cohort setup and sample-list preparation have already been completed. It currently expects versioned input and output paths to be supplied by the user.

The notebook currently performs the following operations:

1. Restricts the normalized proteomics matrix to the AoU sample list.
2. Calculates assay- and participant-level missingness and retains participants with complete assay measurements.
3. Uses robust sample connectivity to identify proteomic outliers.
4. Maps UniProt identifiers to Ensembl identifiers and GENCODE v48 GRCh38 gene coordinates.
5. Rank-normalizes protein measurements within each ancestry group.
6. Calculates ancestry-specific proteomic phenotype PCs.
7. Combines proteomic phenotype PCs with genetic PCs and writes ancestry-specific covariate files.

The notebook's exact file paths, release identifiers, sample counts, and output locations should be captured alongside each analysis run. The notebook is a preparation artifact for the current process; it is not yet a finalized standalone production workflow.

## Shared Genotype Preparation

Proteomics pQTL analyses use the same genotype preparation components as the AoU RNA-seq QTL workflows:

* [`PrepareGenotypes.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/PrepareGenotypes.wdl) filters and exports genotype data, creates PLINK2 files, and calculates genetic PCs.
* [`calculateGenotypeDosage.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/prepare_QTL/calculateGenotypeDosage.wdl) prepares genotype dosage inputs.
* The proteomics covariate notebook combines genetic PCs with proteomic phenotype PCs. The resulting covariate paths are then added to the mapping-input table manually; the shared [`MergeCovariates.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/prepare_QTL/MergeCovariates.wdl) is not listed as part of the current proteomics run.

The genotype files supplied to TensorQTL are the PLINK2 `pgen`, `pvar`, and `psam` outputs from genotype preparation. The sample identifiers in the genotype, phenotype, and covariate files must match before association testing.

## TensorQTL Association

The shared [`tensorqtl_cis_permutations.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/tensorQTL_cis_permutations/tensorqtl_cis_permutations.wdl) workflow runs TensorQTL in cis mode with:

* PLINK2 genotype files
* A proteomics phenotype BED file
* Merged genetic and proteomic covariates
* An analysis prefix
* Optional FDR, p-value, q-value lambda, seed, and command-line flag parameters

The workflow emits a compressed cis-QTL result file and a TensorQTL log. The current WDL also accepts optional phenotype groups for sQTL analyses; phenotype groups are generally not required for the proteomics pQTL use case.

## Fine-Mapping and Aggregation

The downstream fine-mapping process is intended to reuse the RNA-seq QTL components. After TensorQTL results are reviewed and filtered, the current plan is to prepare mapping inputs, run [`susieR_workflow.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/susieR_workflow.wdl) for each ancestry and for the combined cohort, calculate allele frequencies with [`calculateAF.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/prepare_QTL/calculateAF.wdl), and aggregate the annotated SuSiE results with [`AggregateSusieWorkflow.wdl`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/AggregateSusieWorkflow.wdl).

Fine-mapping, aggregation, reporting, and functional-enrichment outputs should be considered in progress until the corresponding workflows and release artifacts are finalized.

## Practical Run Notes

* Run the cohort setup and ancestry partitioning before generating genotype or proteomic phenotype inputs.
* Use the same sample identifiers and population definition for the genotype, phenotype, and covariate files in each analysis.
* Record the proteomics release, ancestry assignment version, GENCODE/Ensembl versions, input paths, output paths, and software/container versions.
* Review missingness, outlier counts, sample overlap, and post-merge missing values before submitting TensorQTL.
* Treat the current workflows and this documentation as in progress until the interfaces and release process are finalized.

## Acknowledgements

The foundational proteomics QTL analyses and supporting scripts were developed by the **Dr. Stephen Montgomery Lab** at Stanford University. Related analysis scripts and workflows are available in the [AoU Multiomics Analysis repository](https://github.com/AoU-Multiomics-Analysis). Special thanks to:

* **Evin Padhi**
* **Jon Nguyen**

and to the other members of the Stanford group who contributed to developing the proteomics processing, covariate preparation, QTL analysis, and fine-mapping methods.

Additional integration, optimization, and workflow migration were performed by the All of Us Multiomics and Broad Pipeline Development teams as part of the WARP workflow suite.

## Feedback

Please help us make these workflows better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions and questions.

</div>
