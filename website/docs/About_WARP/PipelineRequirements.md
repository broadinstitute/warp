---
sidebar_position: 2
---

# Pipeline Requirements

All released WARP pipelines should meet the following criteria:

1. Pipelines are written in WDL, modular, and optimized for Google cloud. The pipelines are often portable to other environments.
2. Pipelines have plumbing (fast) tests and scientific tests for catching any unintended changes in data processing.
3. Pipelines are semantically versioned, have a changelog, and are packaged into a release with all dependencies.
4. Pipelines are released to Dockstore automatically upon release and are available in the cloud-based [Terra platform](https://app.terra.bio/).
5. Pipelines have example inputs alongside the pipeline.
6. Pipelines have a Readme.md describing the pipeline.
7. Pipelines use public docker containers and only open source tools.
8. Pipelines are developed in a collaboration between software engineers and a scientific owner of each pipeline who must approve all changes to the code and verify that any resulting changes in outputs are scientifically valid.
9. Pipelines have a publication-style methods section to enable easy citation.

:::tip
Read more about our pipeline development in our [Best Practices](./BestPractices.md) documentation.
:::

When citing WARP, please use the following:

Degatano, Kylee, Aseel Awdeh, Robert Sidney Cox III, Wes Dingman, George Grant, Farzaneh Khajouei, Elizabeth Kiernan, et al. 2025. "Warp Analysis Research Pipelines: Cloud-Optimized Workflows for Biological Data Processing and Reproducible Analysis." _Bioinformatics (Oxford, England)_, September, https://doi.org/10.1093/bioinformatics/btaf494
