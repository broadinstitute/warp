---
sidebar_position: 1
slug: /Pipelines/Smart-seq2_Multi_Sample_Pipeline/README
---

# Smart-seq2 Multi-Sample Overview
:::warning
9/12/2014

We are deprecating the Smart-seq2 Multi-Sample Pipeline. Although the code will continue to be available, we are no longer supporting it. For an alternative, see the [Smart-seq2 Single Nucleus Multi Sample workflow](../../Pipelines/Smart-seq2_Single_Nucleus_Multi_Sample_Pipeline/README.md).
:::

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [MultiSampleSmartSeq2_v2.2.21](https://github.com/broadinstitute/warp/releases) | December, 2023 | Elizabeth Kiernan | Please [file an issue in WARP](https://github.com/broadinstitute/warp/issues). |

## Introduction

The Smart-seq2 Multi-Sample (Multi-SS2) Pipeline is a wrapper around the [Smart-seq2 Single Sample](../Smart-seq2_Single_Sample_Pipeline/README) pipeline. It is developed by the Data Coordination Platform of the Human Cell Atlas to process single-cell RNAseq (scRNAseq) data generated by Smart-seq2 assays. The workflow processes multiple cells by importing and running the [Smart-seq2 Single Sample workflow](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/smartseq2_single_sample/SmartSeq2SingleSample.wdl) for each cell (sample) and then merging the resulting Loom matrix output into a single Loom matrix containing raw counts and TPMs.

Full details about the Smart-seq2 Pipeline can be read in the [Smart-seq2 Single Sample Overview](https://broadinstitute.github.io/warp/docs/Pipelines/Smart-seq2_Single_Sample_Pipeline/README) in GitHub.

The Multi-SS2 workflow can also be run in [Terra](https://app.terra.bio), a cloud-based analysis platform. The Terra [Smart-seq2 public workspace](https://app.terra.bio/#workspaces/featured-workspaces-hca/HCA%20Smart-seq2%20Multi%20Sample%20Pipeline) contains the Smart-seq2 workflow, workflow configurations, required reference data and other inputs, and example testing data.

:::tip Want to use the Multi-SS2 Pipeline for your publication?
Check out the [Smart-seq2 Publication Methods](./smart-seq2.methods.md) to get started!
:::

## Inputs

There are two example configuration (JSON) files available for testing the Multi-SS2 workflow. Both examples are also preloaded in the Terra [Smart-seq2 public workspace](https://app.terra.bio/#workspaces/featured-workspaces-hca/HCA%20Smart-seq2%20Multi%20Sample%20Pipeline).
* [human_single_example.json](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/smartseq2_multisample/human_single_example.json): Configurations for an example single-end human dataset consisting of two samples (cells)
* [mouse_paired_example.json](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/smartseq2_multisample/mouse_paired_example.json): Configurations for an example paired-end mouse dataset consisting of two samples (cells)


### Sample and Reference Inputs

The workflow’s sample inputs are listed in the table below. Reference inputs are identical to those
specified in the [Smart-seq2 Single Sample Overview](../Smart-seq2_Single_Sample_Pipeline/README).

The workflow processes both single- and paired-end samples; however, these samples can not be mixed in the same run.

| Input name | Input Description | Input Type |
| --- | --- | --- |
| fastq1_input_files | Cloud locations for each read1 file | Array of strings |
| fastq2_input_files | Optional cloud locations for each read2 file if running paired-end samples |Array of strings |
| input_ids | Unique identifiers or names for each cell; can be a UUID or human-readable name | Array of strings |
| input_names | Optional unique identifiers/names to further describe each cell. If `input_id` is a UUID, the `input_name` could be used as a human-readable identifier | String |
| batch_id | Identifier for the batch of multiple samples | String |
| batch_name | Optional string to describe the batch or biological sample | String |
| input_name_metadata_field | Optional input describing, when applicable, the metadata field containing the `input_names` | String |
| input_id_metadata_field | Optional string describing, when applicable, the metadata field containing the `input_ids` | String |
| `project_id` | Optional project identifier; usually a number | String |
| `project_name` | Optional project identifier; usually a human-readable name | String |
| `library` | Optional description of the sequencing method or approach | String |
| `organ` | Optional description of the organ from which the cells were derived | String |
| `species` | Optional description of the species from which the cells were derived | String |
| `paired-end` | Boolean for whether samples are paired-end or not | Boolean |


### Additional Input

The reference inputs are identical to those specified in the "Additional Reference Inputs" section of the [Smart-seq2 Single Sample Overview](https://broadinstitute.github.io/warp/docs/Pipelines/Smart-seq2_Single_Sample_Pipeline/README#inputs).


### Smart-seq2 Multi-Sample Task Summary

The Multi-SS2 Pipeline calls two tasks:

1) [SmartSeq2SingleSample](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/smartseq2_single_sample/SmartSeq2SingleSample.wdl): a task that runs the Smart-seq2 Single Sample workflow
2) [SmartSeq2PlateAggregation](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/LoomUtils.wdl): the wrapper pipeline that aggregates the results


### Outputs

| Output file name | Output Description | Output Type |
| --- | --- | --- |
| bam_files | An array of genome-aligned BAM files (one for each sample) generated with HISAT2  | Array |
| bam_index_files |  An array of BAM index files generated with HISAT2 | Array |
| loom_output | A single Loom cell-by-gene matrix containing raw counts and TPMs for every cell  | File |

The final Loom matrix is an aggregate of all the individual Loom matrices generated using the [Smart-seq2 Single Sample workflow](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/smartseq2_single_sample/SmartSeq2SingleSample.wdl).

The aggregated Loom filename contains the `batch_id` prefix, which is the string specified in the input configuration. The `batch_id` is also set as a global attribute in the Loom.

Both the individual sample Loom files and individual BAM files are described in the [Smart-seq2 Single Sample Overview](../Smart-seq2_Single_Sample_Pipeline/README).

:::warning Zarr Array Deprecation Notice June 2020
Please note that we have deprecated the previously used Zarr array output. The pipeline now uses the Loom file format as the default output.
:::

## Validation
The Multi-SS2 Pipeline has been validated for processing human and mouse, stranded or unstranded, paired- or single-end, and plate- or fluidigm-based Smart-seq2 datasets (see links to validation reports in the table below).

| Workflow Configuration | Link to Report |
| --- | --- |
| Mouse paired-end | [Report](https://docs.google.com/document/d/12zGTFROrcXEByt9z0h06qjSqb9vWutn28Tx6YiND1Ds/edit)
| Human and mouse single-end | [Report](https://docs.google.com/document/d/1MonsTG8UnROHZ_XpulrSZNTxO988KEH6T6h45plFYQg/edit#heading=h.ixoqmhbabdvh) |
| Human stranded fluidigm | [Report](https://docs.google.com/document/d/1FEg86Tlu657j9Kjw_v3keFQRXcBIs8gOqCwLbPSP-C0/edit#heading=h.wjr8otl7zg14) |

## Versioning

Release information for the Multi-SS2 Pipeline can be found in the [changelog](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/smartseq2_multisample/MultiSampleSmartSeq2.changelog.md). Please note that any major changes to the Smart-seq2 pipeline will be documented in the [Smart-seq2 Single Sample changelog](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/smartseq2_single_sample/SmartSeq2SingleSample.changelog.md).

## Citing the Smart-seq2 Multi-Sample Pipeline

If you use the Smart-seq2 Multi-Sample Pipeline in your research, please identify the pipeline in your methods section using the [Smart-seq2 Multi-Sample SciCrunch resource identifier](https://scicrunch.org/resources/data/record/nlx_144509-1/SCR_018920/resolver?q=SCR_018920&l=SCR_018920&i=rrid:scr_018920).

* Ex: *Smart-seq2 Multi-Sample Pipeline (RRID:SCR_018920)*

Please also consider citing our preprint:

Degatano, K.; Awdeh, A.; Dingman, W.; Grant, G.; Khajouei, F.; Kiernan, E.; Konwar, K.; Mathews, K.; Palis, K.; Petrillo, N.; Van der Auwera, G.; Wang, C.; Way, J.; Pipelines, W. WDL Analysis Research Pipelines: Cloud-Optimized Workflows for Biological Data Processing and Reproducible Analysis. Preprints 2024, 2024012131. https://doi.org/10.20944/preprints202401.2131.v1

## Consortia Support
This pipeline is supported and used by the [Human Cell Atlas](https://www.humancellatlas.org/) (HCA) project. 

If your organization also uses this pipeline, we would love to list you! Please reach out to us by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues).

## Have Suggestions?
Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

