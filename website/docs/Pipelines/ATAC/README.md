---
sidebar_position: 1
---

# ATAC Overview

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [not release](https://github.com/broadinstitute/warp/releases) | May, 2023 | Kaylee Mathews | Please file GitHub issues in warp or contact [the WARP team](mailto:warp-pipelines-help@broadinstitute.org) |

![Multiome_diagram]()

## Introduction to the ATAC workflow
ATAC is an open-source, cloud-optimized pipeline developed in collaboration with members of the [BRAIN Initiative Cell Census Network](https://biccn.org/) (BICCN), including the Allen Institute, NeMO, Kai Zhang (SnapATAC2), and Alex Dobin. It supports the processing of 10x single-nucleus data generated with ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing), a technique used in molecular biology to assess genome-wide chromatin accessibility. 

generated with the [10x Genomics Multiome assay](https://www.10xgenomics.com/products/single-cell-multiome-atac-plus-gene-expression).

The pipeline corrects cell barcodes, aligns reads to the genome, and producesa fragment file as well as per barcode metrics. 


## Quickstart table
The following table provides a quick glance at the Multiome pipeline features:

| Pipeline features | Description | Source |
|--- | --- | --- |
| Assay type | 10x single cell or single nucleus gene expression (GEX) and ATAC | [10x Genomics](https://www.10xgenomics.com)
| Overall workflow  | Barcode correction, read alignment, and fragment quanitification |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Genomic Reference Sequence | GRCh38 human genome primary sequence | GENCODE |
| Aligner | BWA-mem | [Dobin, et al.,2021](https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1) |
| Fragment quantification | SnapATAC2 |
| Data input file format | File format in which sequencing data is provided | [FASTQ](https://academic.oup.com/nar/article/38/6/1767/3112533) |
| Data output file format | File formats in which Multiome output is provided | BAM, h5ad |

## Set-up

### ATAC installation

To download the latest ATAC release, see the release tags prefixed with "Multiome" on the WARP [releases page](https://github.com/broadinstitute/warp/releases). All Multiome pipeline releases are documented in the [Multiome changelog](). 

To discover and search releases, use the WARP command-line tool [Wreleaser](https://github.com/broadinstitute/warp/tree/master/wreleaser).


ATAC can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH compliant, flexible workflow management system that supports multiple computing platforms. The workflow can also be run in [Terra](https://app.terra.bio), a cloud-based analysis platform. 



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
| whitelist                 | Whitelist file for cellular barcodes                                                                 |
| adapter_seq_read1         | TrimAdapters input: Sequence adapter for read 1 fastq                                                |
| adapter_seq_read3         | TrimAdapters input: Sequence adapter for read 3 fastq                                                |

## ATAC tasks and tools

Overall, the Optimus workflow:
1. Corrects CBs and partitions FASTQs by CB.
1. Aligns reads.
1. Merges aligned BAMs
1. Generates a fragment file
1. Calculates per cell barcode fragment metrics.


The tools each ATAC task employs are detailed in the table below. 

To see specific tool parameters, select the task WDL link in the table; then view the `command {}` section of the task in the WDL script. To view or use the exact tool software, see the task's Docker image which is specified in the task WDL `# runtime values` section as `String docker =`.

| Task name and WDL link | Tool | Software | Description | 
| --- | --- | --- | ------------------------------------ | 
| [FastqProcessing as SplitFastq]() | fastqprocess | custom | Corrects cell barcodes and splits fastqs into 10 GB FASTQs that are grouped by cell barcode with each read having the corrected and raw barcode in the read name. |
| [TrimAdapters]() | Cutadapt | cutadapt | Trims adaptor sequences. |
| [BWAPairedEndAlignment]() | bwa-mem | bwa-mem | Aligns reads from each set of partitioned FASTQ files to the genome and outputs a BAM. |
| [AddCBtags]() | view, sort | samtools | Partitions the barcodes from the BAM readname into the CB (corrected barcodes) and CR tags (raw barcodes). |
| [Merge.MergeSortBamFiles as MergeBam]() | MergeSamFiles | Picard | Merges each BAM into a final aligned BAM with corrected cell barcodes in the CB tag. |
| [CreateFragmentFile]() | make_fragment_file, import_data | SnapATAC2 | Generates a fragment file from the final aligned BAM and outputs quality metrics in h5ad. |


## Output Variables

| Variable Name         | Description                                    |
|----------------------|------------------------------------------------|
| bam_aligned_output    | Aligned BAM file                               |
| fragment_file         | Fragment file                                  |
| snap_metrics          | Snap metrics file                              |