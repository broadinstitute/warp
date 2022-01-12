---
sidebar_position: 1
---

# RNA with UMIs Overview

| Pipeline Version | Date Updated | Documentation Authors | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [RNAWithUMIsPipeline_v0.1.0](https://github.com/broadinstitute/warp/releases?q=RNAwithUMIs&expanded=true) | January, 2022 | [Elizabeth Kiernan](mailto:ekiernan@broadinstitute.org) & [Kaylee Mathews](mailto:kmathews@broadinstitute.org)| Please file GitHub issues in warp or contact [Kylee Degatano](mailto:kdegatano@broadinstitute.org) |

<!--- Pipeline diagram will go here --->

## Introduction to the RNA with UMIs workflow

The [RNA with UMIs Pipeline](https://github.com/broadinstitute/warp/blob/develop/pipelines/broad/rna_seq/RNAWithUMIsPipeline.wdl) is an open-source, cloud-optimized workflow for processing total RNA isolated with the Transcriptome Capture (TCap) method. TCap is a technique that hybridizes exome baits to cDNA library preparations to better facilitate RNA-sequencing in low-input or poor quality (degraded) samples. These libraries may additionally be prepared with Unique Molecular Identifiers (UMIs) which can help distinguish biological signal from noise resulting from PCR amplification. 

Overall, the workflow performs UMI correction, aligns reads to the genome, quantifies gene counts, and calculates quality metrics. The workflow produces genome- and transcriptome-aligned BAMs with indices and a merged quality metrics file. 

While this workflow was created to be used with TCap RNA-seq data, it can be used to process any bulk RNA-seq data. 

<!--- tip for methods section will go here --->

<!--- quickstart table will go here --->

## Set-up

### RNA with UMIs pipeline installation

To download the latest release of the RNA with UMIs pipeline, see the release tags prefixed with "RNAwithUMIs" on the WARP [releases page](https://github.com/broadinstitute/warp/releases). All releases of the RNA with UMIs pipeline are documented in the [RNA with UMIs changelog](https://github.com/broadinstitute/warp/blob/develop/pipelines/broad/rna_seq/RNAWithUMIsPipeline.changelog.md). 

To search releases of this and other pipelines, use the WARP command-line tool [Wreleaser](https://github.com/broadinstitute/warp/tree/develop/wreleaser).

The RNA with UMIs pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH compliant, flexible workflow management system that supports multiple computing platforms. The workflow can also be run in [Terra](https://app.terra.bio), a cloud-based analysis platform. <!--- link to public workspace will go here --->