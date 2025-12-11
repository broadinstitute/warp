---
sidebar_position: 2
---

# Smart-seq2 Multi-Sample v2.2.21 Publication Methods

Below we provide an example methods section for a publication. For the complete pipeline documentation, see the [Smart-seq2 Multi-Sample Overview](./README.md).

## Methods

Data preprocessing and count matrix construction for a sample batch (or plate) were performed using the Smart-seq2 Multi-Sample v2.2.21 Pipeline (RRID:SCR_018920). For each cell in the batch, paired- or single-end FASTQ files were first processed with the Smart-seq2 Single Sample v5.1.20 Pipeline (RRID:SCR_021228). Reads were aligned to the GENCODE mouse (M21) or human (V27) reference genome using HISAT2 v2.1.0 with default parameters in addition to `--k 10` options. Metrics were collected and duplicate reads marked using the Picard v.2.26.10 `CollectMultipleMetrics` and `CollectRnaSeqMetrics`, and MarkDuplicates functions with validation_stringency=silent. For transcriptome quantification, reads were aligned to the GENCODE transcriptome using HISAT2 v2.1.0 with `--k 10 --no-mixed  --no-softclip  --no-discordant --rdg 99999999,99999999 --rfg 99999999,99999999 --no-spliced-alignment` options. Gene expression was calculated using RSEM v1.3.0’s `rsem-calculate-expression --calc-pme --single-cell-prior`. QC metrics, RSEM TPMs and RSEM estimated counts were exported to a single Loom file for each cell. All individual Loom files for the entire batch were aggregated into a single Loom file for downstream processing. The final output included the unfiltered Loom and the tagged, unfiltered individual BAM files.

An example of the pipeline and outputs can be found in [Terra](https://app.terra.bio/#workspaces/featured-workspaces-hca/HCA%20Smart-seq2%20Multi%20Sample%20Pipeline) and additional documentation can be found in the [Smart-seq2 Multi-Sample Overview](./README.md). Examples of genomic references, whitelists, and other inputs are available in the warp repository (see the *_example.json files at [human_single_example](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/smartseq2_multisample/human_single_example.json).
