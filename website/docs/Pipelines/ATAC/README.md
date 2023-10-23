---
sidebar_position: 1
slug: /Pipelines/ATAC/README

---

# ATAC Overview

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [1.1.1](https://github.com/broadinstitute/warp/releases) | October, 2023 | Kaylee Mathews | Please file GitHub issues in warp or contact [the WARP team](mailto:warp-pipelines-help@broadinstitute.org) |

## Introduction to the ATAC workflow
ATAC is an open-source, cloud-optimized pipeline developed collaboration with members of the [BRAIN Initiative](https://braininitiative.nih.gov/) (BICCN and [BICAN](https://brainblog.nih.gov/brain-blog/brain-issues-suite-funding-opportunities-advance-brain-cell-atlases-through-centers) Sequencing Working Group) and [SCORCH](https://nida.nih.gov/about-nida/organization/divisions/division-neuroscience-behavior-dnb/basic-research-hiv-substance-use-disorder/scorch-program) (see [Acknowledgements](#acknowledgements) below). It supports the processing of 10x single-nucleus data generated with 10x Multiome [ATAC-seq (Assay for Transposase-Accessible Chromatin)](https://www.10xgenomics.com/products/single-cell-multiome-atac-plus-gene-expression), a technique used in molecular biology to assess genome-wide chromatin accessibility. 

This workflow is the ATAC component of the [Mutiome wrapper workflow](../Multiome_Pipeline/README). It corrects cell barcodes (CBs), aligns reads to the genome, and produces a fragment file as well as per barcode metrics. 


## Quickstart table
The following table provides a quick glance at the ATAC pipeline features:

| Pipeline features | Description | Source |
|--- | --- | --- |
| Assay type | 10x single cell or single nucleus ATAC | [10x Genomics](https://www.10xgenomics.com)
| Overall workflow  | Barcode correction, read alignment, and fragment quanitification |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Genomic Reference Sequence | GRCh38 human genome primary sequence | GENCODE |
| Aligner | BWA-mem | [Li H. and Durbin R., 2009](http://www.ncbi.nlm.nih.gov/pubmed/19451168) |
| Fragment quantification | SnapATAC2 | [Zhang, K. et al., 2021](https://pubmed.ncbi.nlm.nih.gov/34774128/)
| Data input file format | File format in which sequencing data is provided | [FASTQ](https://academic.oup.com/nar/article/38/6/1767/3112533) |
| Data output file format | File formats in which ATAC output is provided | TSV, h5ad, BAM |

## Set-up

### ATAC installation

To download the latest ATAC release, see the release tags prefixed with "Multiome" on the WARP [releases page](https://github.com/broadinstitute/warp/releases). All ATAC pipeline releases are documented in the [ATAC changelog](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/multiome/atac.changelog.md). 

To discover and search releases, use the WARP command-line tool [Wreleaser](https://github.com/broadinstitute/warp/tree/master/wreleaser).

ATAC can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH compliant, flexible workflow management system that supports multiple computing platforms. The workflow can also be run in [Terra](https://app.terra.bio), a cloud-based analysis platform. 

## Input Variables
The following describes the inputs of the ATAC workflow. For more details on how default inputs are set for the Multiome workflow, see the [Multiome overview](../Multiome_Pipeline/README).

| Variable name | Description |
| --- | --- |
| read1_fastq_gzipped | Fastq inputs (array of compressed read 1 FASTQ files). |
| read2_fastq_gzipped | Fastq inputs (array of compressed read 2 FASTQ files containing cellular barcodes). |
| read3_fastq_gzipped | Fastq inputs (array of compressed read 3 FASTQ files). |
| output_base_name | Output prefix/base name for all intermediate files and pipeline outputs. |
| tar_bwa_reference | BWA reference (tar file containing reference fasta and corresponding files). |
| atac_gtf | CreateFragmentFile input variable: GTF file for SnapATAC2 to calculate TSS sites of fragment file.|
| chrom_sizes | CreateFragmentFile input variable: Text file containing chrom_sizes for genome build (i.e., hg38) |
| whitelist | Whitelist file for ATAC cellular barcodes. |
| adapter_seq_read1 | TrimAdapters input: Sequence adapter for read 1 fastq. |
| adapter_seq_read3 | TrimAdapters input: Sequence adapter for read 3 fastq. |

## ATAC tasks and tools

Overall, the ATAC workflow:
1. Corrects CBs and partitions FASTQs by CB.
1. Aligns reads.
1. Merges aligned BAMs
1. Generates a fragment file
1. Calculates per cell barcode fragment metrics.

The tools each ATAC task employs are detailed in the table below. 

To see specific tool parameters, select the task WDL link in the table; then view the `command {}` section of the task in the WDL script. To view or use the exact tool software, see the task's Docker image which is specified in the task WDL `# runtime values` section as `String docker =`.

| Task name and WDL link | Tool | Software | Description | 
| --- | --- | --- | ------------------------------------ | 
| [FastqProcessing as SplitFastq](https://github.com/broadinstitute/warp/blob/develop/tasks/skylab/FastqProcessing.wdl) | fastqprocess | custom | Dynamically selects the correct barcode orientation, corrects cell barcodes and splits FASTQs into 30 GB FASTQs that are grouped by cell barcode with each read having the corrected (CB) and raw barcode (CR) in the read name. |
| [TrimAdapters](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/multiome/atac.wdl) | Cutadapt v4.4 | cutadapt | Trims adaptor sequences. |
| [BWAPairedEndAlignment](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/multiome/atac.wdl) | bwa-mem2 | mem | Aligns reads from each set of partitioned FASTQ files to the genome and outputs a BAM with ATAC barcodes in the CB:Z tag. |
| [Merge.MergeSortBamFiles as MergeBam](https://github.com/broadinstitute/warp/blob/develop/tasks/skylab/MergeSortBam.wdl) | MergeSamFiles | Picard | Merges each BAM into a final aligned BAM with corrected cell barcodes in the CB tag. |
| [CreateFragmentFile](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/multiome/atac.wdl) | make_fragment_file, import_data | SnapATAC2 | Generates a fragment file from the final aligned BAM and outputs per barcode quality metrics in h5ad. |


## Output variables

| Output variable name | Filename, if applicable | Output format and description |
|--- | --- | --- | 
| bam_aligned_output | `<input_id>`.bam | BAM containing aligned reads from ATAC workflow. |
| fragment_file | `<input_id>`.fragments.tsv | TSV containing fragment start and stop coordinates per barcode. | 
| snap_metrics | `<input_id`.metrics.h5ad | h5ad (Anndata) containing per barcode metrics from SnapATAC2. |

## Acknowledgements

We are immensely grateful to the members of the BRAIN Initiative (BICAN Sequencing Working Group) and SCORCH for their invaluable and exceptional contributions to this pipeline. Our heartfelt appreciation goes to Alex Dobin, Aparna Bhaduri, Alec Wysoker, Anish Chakka, Brian Herb, Daofeng Li, Fenna Krienen, Guo-Long Zuo, Jeff Goldy, Kai Zhang, Khalid Shakir, Bo Li, Mariano Gabitto, Michael DeBerardine, Mengyi Song, Melissa Goldman, Nelson Johansen, James Nemesh, and Theresa Hodges for their unwavering dedication and remarkable efforts.  