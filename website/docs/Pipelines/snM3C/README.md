---
sidebar_position: 1
slug: /Pipelines/snM3C/README
---
# Single Nucleus Methyl-Seq and Chromatin Capture (snM3C) Overview

| Pipeline Version | Date Updated | Documentation Authors | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [snM3C_v1.0.0](https://github.com/broadinstitute/warp/releases) | August, 2023 | [Kaylee Mathews](mailto:warp-pipelines-help@broadinsitute.org) | Please file GitHub issues in the [WARP repository](https://github.com/broadinstitute/warp/issues) |


## Introduction to snM3C

The Single Nucleus Methly-Seq and Chromatin Capture (snM3C) workflow is a cloud-based computational workflow for processing single-nucleus methylome and chromatin contact (snM3C) sequencing data. The workflow is designed to demultiplex raw sequencing reads, align them, call chromatin contacts, and generate summary metrics. It is developed in collaboration Hanqing Liu and the laboratory of Joseph Ecker. For more information about the snM3C tools and analysis, please see the [YAP documentation](https://hq-1.gitbook.io/mc/) or the [cemba_data](https://github.com/lhqing/cemba_data) GitHub repository created by Hanqing Liu.

## Set-up

### Installation

To use the latest release of the snM3C pipeline, visit the [WARP releases page](https://github.com/broadinstitute/warp/releases) and download the desired version.

<!--- Add a comment about running an old version of the workflow --->

### Running the Workflow

To download the latest release of the snM3C pipeline, see the release tags prefixed with "snM3C" on the WARP [releases page](https://github.com/broadinstitute/warp/releases). All releases of the snM3C pipeline are documented in the [snM3C changelog](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/snM3C/snM3C.changelog.md). 

To search releases of this and other pipelines, use the WARP command-line tool [Wreleaser](https://github.com/broadinstitute/warp/tree/develop/wreleaser).

<!--- add a comment about running an old version of the workflow --->

The snM3C pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH compliant, flexible workflow management system that supports multiple computing platforms. The workflow can also be run in [Terra](https://app.terra.bio), a cloud-based analysis platform. 

### Inputs

The snM3C workflow requires a JSON configuration file specifying the input files and parameters for the analysis. An example configuration file can be found in the [snM3C direcotry](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/snM3C/snM3C_inputs.json).

The main input files and parameters include:

| Parameter | Description |
| ---| --- |
| fastq_input_read1 | Array of multiplexed FASTQ files for read 1 |
| fastq_input_read2 | Array of multiplexed FASTQ files for read 2 |
| random_primer_indexes | File containing random primer indexes |
| plate_id | String specifying the plate ID |
| output_basename | String specifying a basename to be used for naming files |
| tarred_index_files | File containing tarred index files for hisat-3 mapping |
| mapping_yaml | File containing YAML configuration for mapping steps with snakemake | 
| snakefile | File containing the snakefile for mapping |
| chromosome_sizes | File containing chromosome sizes information |
| genome_fa | File containing the reference genome in FASTA format | 


## Tasks and Tools
The workflow contains two tasks described below. The parameters and more details about these tools can be found in the [YAP documentation](https://hq-1.gitbook.io/mc/).

| Task name | Tool | Software | Description |
| --- | --- | --- | --- |
| Demultiplexing | cutadapt | cutadapt | Performs demultiplexing to cell-level FASTQ files |
| Mapping | hisat-3 | hisat-3 | Performs trimming, alignment and calling chromatin contacts with a [custom snakemake](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/snM3C/Config%20files/Snakemake-file/Snakefile) file developed by Hanqing Liu. |

## Outputs

The snM3C workflow produces the following main outputs:

| Output | Description | 
| ---| --- |
| mappingSummary | Mapping summary file in CSV format |
| allcFiles | Tarred file containing allc files |
| allc_CGNFiles| Tarred file containing CGN context-specific allc files | 
| bamFiles | Tarred file containing cell-level aligned BAM files |
| detail_statsFiles | Tarred file containing detail stats files | 
| hicFiles | Tarred file containing Hi-C files |


## Versioning

All snM3C pipeline releases are documented in the [pipeline changelog](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/snM3C/snM3C.changelog.md).

<!--- Citing the pipeline will go here --->

## Feedback

For questions, suggestions, or feedback related to the snM3C pipeline, please contact [the WARP team](mailto:warp-pipelines-help@broadinstitute.org). Your feedback is valuable for improving the pipeline and addressing any issues that may arise during its usage.

<!--- Validation will go here --->



<!--- FAQs will go here --->


