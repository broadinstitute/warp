---
sidebar_position: 1
---

# Multiome Overview

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [not release](https://github.com/broadinstitute/warp/releases) | May, 2023 | Elizabeth Kiernan | Please file GitHub issues in warp or contact [the WARP team](mailto:warp-pipelines-help@broadinstitute.org) |

![Multiome_diagram]()

## Introduction to the Multiome workflow
Multiome is an open-source, cloud-optimized pipeline developed in collaboration with members of the [BRAIN Initiative Cell Census Network](https://biccn.org/) (BICCN), including the Allen Institute, NeMO, Kai Zhang (SnapATAC2), and Alex Dobin. It supports the processing of 10x 3' single-cell and single-nucleus count data generated with the [10x Genomics Multiome assay](https://www.10xgenomics.com/products/single-cell-multiome-atac-plus-gene-expression).

The pipeline consists of two subworkflows: the Optimus workflow for GEX data and the ATAC workflow for single-cell ATAC data. The GEX component performs barcode and UMI correction, aligns reads to the genome, and produces both quality metrics per barcode and gene and a raw cell-by-gene count matrix. The ATAC component corrects cell barcodes, aligns reads to the genome, and producesa fragment file as well as per barcode metrics. 

Read more about the GEX workflow in the Optimus overview for the ATAC workflow in the ATAC workflow. 


## Quickstart table
The following table provides a quick glance at the Multiome pipeline features:

| Pipeline features | Description | Source |
|--- | --- | --- |
| Assay type | 10x single cell or single nucleus gene expression (GEX) and ATAC | [10x Genomics](https://www.10xgenomics.com)
| Overall workflow  | Barcode correction, read alignment, gene and fragment quanitification |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Genomic Reference Sequence | GRCh38 human genome primary sequence | GENCODE |
| Gene annotation reference (GTF) | Reference containing gene annotations | Gencode |
| Aligners | STARsolo (GEX), BWA-mem (ATAC) | [Dobin, et al.,2021](https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1) |
| Transcript and fragment quantification | STARsolo (GEX), SnapATAC2 (ATAC)
| Data input file format | File format in which sequencing data is provided | [FASTQ](https://academic.oup.com/nar/article/38/6/1767/3112533) |
| Data output file format | File formats in which Multiome output is provided | BAM, h5ad |

## Set-up

### Multiome installation

To download the latest Multiome release, see the release tags prefixed with "Multiome" on the WARP [releases page](https://github.com/broadinstitute/warp/releases). All Multiome pipeline releases are documented in the [Multiome changelog](). 

To discover and search releases, use the WARP command-line tool [Wreleaser](https://github.com/broadinstitute/warp/tree/master/wreleaser).

If you’re running an ultiomeus workflow version prior to the latest release, the accompanying documentation for that release may be downloaded with the source code on the WARP [releases page](https://github.com/broadinstitute/warp/releases) (see the source code folder).

Multiome can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH compliant, flexible workflow management system that supports multiple computing platforms. The workflow can also be run in [Terra](https://app.terra.bio), a cloud-based analysis platform. 


### Single Cell ATAC Inputs





# Workflow: ATAC

**Version:** 1.0

**Description:** Processing for single-cell ATAC-seq data from the level of raw fastq reads. This is the first step of the multiome pipeline. ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) is a technique used in molecular biology to assess genome-wide chromatin accessibility. This pipeline processes 10x Genomics Multiome ATAC FASTQ files.

## Input Variables

| Variable Name              | Description                                                                                          |
|---------------------------|------------------------------------------------------------------------------------------------------|
| read1_fastq_gzipped       | Fastq inputs (array of compressed read 1 FASTQ files)                                                |
| read2_fastq_gzipped       | Fastq inputs (array of compressed read 2 FASTQ files containing cellular barcodes)                   |
| read3_fastq_gzipped       | Fastq inputs (array of compressed read 3 FASTQ files)                                                |
| output_base_name          | Output prefix/base name for all intermediate files and pipeline outputs                              |
| tar_bwa_reference         | BWA reference (tar file containing reference fasta and corresponding files)                          |
| barcodes_in_read_name     | CreateFragmentFile input variable: Boolean indicating whether barcodes are in read names             |
| atac_gtf                  | CreateFragmentFile input variable: GTF file for SnapATAC2 to calculate TSS sites of fragment file    |
| chrom_sizes               | CreateFragmentFile input variable: Text file containing chrom_sizes for genome build (i.e., hg38)    |
| monitoring_script         | Script for monitoring tasks                                                                           |
| whitelist                 | Whitelist file for cellular barcodes                                                                 |
| adapter_seq_read1         | TrimAdapters input: Sequence adapter for read 1 fastq                                                |
| adapter_seq_read3         | TrimAdapters input: Sequence adapter for read 3 fastq                                                |

## Output Variables

| Variable Name         | Description                                    |
|----------------------|------------------------------------------------|
| bam_aligned_output    | Aligned BAM file                               |
| fragment_file         | Fragment file                                  |
| snap_metrics          | Snap metrics file                              |