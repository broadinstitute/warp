---
sidebar_position: 1
slug: /Pipelines/SlideTags_Pipeline/README
---

# Slide-tags Pipeline Overview

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| v1.0.0 | March, 2025 | WARP Pipelines | Please [file an issue in WARP](https://github.com/broadinstitute/warp/issues) |

![Slide-tags_diagram](Slide-tags_diagram.png)

## Introduction to the Slide-tags Pipeline

The **Slide-tags Pipeline** is an open-source, cloud-optimized workflow for processing spatial transcriptomics data. It supports data derived from spatially barcoded sequencing technologies, including Slide-tags-based single-molecule profiling. The pipeline processes raw sequencing data into spatially resolved gene expression matrices, ensuring accurate alignment, spatial positioning, and quantification.

This workflow integrates multiple processing steps, including barcode extraction, spatial alignment, transcript counting, and output generation in formats compatible with community tools.

## Quickstart Table

| Pipeline Features | Description | Source |
|--- | --- | --- |
| Assay type | Spatial transcriptomics using Slide-tags | [Macosko Lab](https://macoskolab.com/) |
| Overall workflow  | Barcode extraction, spatial positioning, transcript quantification | Original code available from [GitHub](https://github.com/MacoskoLab/Macosko-Pipelines); WDL workflow available in WARP. |
| Workflow language | WDL | [openWDL](https://github.com/openwdl/wdl) |
| Data input file format | File format in which sequencing data is provided | [FASTQ](https://academic.oup.com/nar/article/38/6/1767/3112533) |
| Data output file format | Output formats for downstream analysis | [HDF5](https://www.hdfgroup.org/) (h5ad), CSV |

## Set-up

### Installation

To download the latest Slide-tags release, see the release tags prefixed with "Optimus" on the WARP [releases page](https://github.com/broadinstitute/warp/releases). All Optimus pipeline releases are documented in the [Slide-tags changelog](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/slidetags/SlideTags.changelog.md).

The pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant workflow manager. Additionally, it can be run in cloud-based analysis platforms such as [Terra](https://app.terra.bio).

### Inputs

The pipeline requires JSON-formatted configuration files detailing input parameters. Required inputs include:

- **Raw FASTQ files** from sequencing
- **Reference genome** and transcript annotation files
- **Spatial barcode whitelist**
- **Spatial positioning reference**

Example input configurations can be found in the `test_inputs` folder of the GitHub repository.

## Slide-tags Pipeline Tasks and Tools

The workflow is composed of several key steps, implemented in separate WDL tasks:

| Task | Tool | Description |
| --- | --- | --- |
| [Optimus](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/Optimus.wdl) | STARsolo | Gene quantification subworkflow that aligns reads to a reference genome and produces a count matrix. Read more in the [Optimus Overview](../Optimus_Pipeline/README.md)| 
| [spatial_count]() | Custom Julia script | Extracts spatial barcodes, performs barcode sequencing error correction, maps reads to spatial barcodes and stores unique (cell, UMI, barcode) triplets in a count matrix, and calculates quality control metrics. Produces an h5 output. |
| [positioning]() | Custom R scripts | Extracts cell barcodes, calculates log-transformed UMI counts, and determines mitochondrial gene percentages. Performs data normalization, PCA, clustering, and UMAP embedding for visualization and produces quality metrics and graphs. Assigns cell barcodes to spatial barcode coordinates. |


Each of these tasks utilizes scripts from the [Macosko Lab Pipelines](https://github.com/MacoskoLab/Macosko-Pipelines) repository, modified for streamlined output handling and stored in [warp-tools]().

## Outputs

| Output | Description | Format |
| ------ | ------ | ------ |
| Aligned BAM | Read alignments with spatial barcodes | BAM |
| Gene Expression Matrix | Cell-by-gene count matrix | h5ad |
| Spatial Barcode Map | Coordinates of detected barcodes | CSV |
| Summary Metrics | Quality control statistics | CSV |

## Versioning

All releases of the pipeline are documented in the repository’s [changelog](https://github.com/MacoskoLab/Macosko-Pipelines/blob/main/CHANGELOG.md).

## Citing the Slide-tags Pipeline

If you use the Slide-tags Pipeline in your research, please cite the original sources:

- **Macosko Lab Pipelines:** https://github.com/MacoskoLab/Macosko-Pipelines

Please also consider citing our preprint:

Degatano, K.; Awdeh, A.; Dingman, W.; Grant, G.; Khajouei, F.; Kiernan, E.; Konwar, K.; Mathews, K.; Palis, K.; Petrillo, N.; Van der Auwera, G.; Wang, C.; Way, J.; Pipelines, W. WDL Analysis Research Pipelines: Cloud-Optimized Workflows for Biological Data Processing and Reproducible Analysis. Preprints 2024, 2024012131. https://doi.org/10.20944/preprints202401.2131.v1

## Acknowledgements

We are immensely grateful Matthew Shabet and the Macosko Lab for development of these analsyes, for their generous time making these scripts FAIR, and for the many hours working with the WARP team to incoporate the scripts into WDL. 

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

