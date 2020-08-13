| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [scATAC 1.1.0 ](scATAC.wdl) | May 18th 2020 | [Elizabeth Kiernan](mailto:ekiernan@broadinstitute.org) | Please file GitHub issues in skylab or contact [Kylee Degatano](mailto:kdegatano@broadinstitute.org) |

- [Overview](#overview)
- [Introduction](#introduction)
  * [Quick Start Table](#quick-start-table)
- [Set-up](#set-up)
  * [Workflow Installation and Requirements](#workflow-installation-and-requirements)
  * [Pipeline Inputs](#pipeline-inputs)
  * [Input File Preparation](#input-file-preparation)
    + [R1 and R2 FASTQ Preparation](#r1-and-r2-fastq-preparation)
    + [Input_reference Preparation](#input_reference-preparation)
- [Workflow Tasks and Tools](#workflow-tasks-and-tools)
  * [Task Summary](#task-summary)
    + [AlignPairedEnd](#alignpairedend)
    + [SnapPre](#snappre)
      - [Filtering Parameters](#filtering-parameters)
      - [Snap QC Metrics](#snap-qc-metrics)
    + [SnapCellByBin](#snapcellbybin)
    + [MakeCompliantBAM](#makecompliantbam)
    + [BreakoutSnap](#breakoutsnap)
- [Outputs](#outputs)
- [Running on Terra](#running-on-terra)
- [Versioning](#versioning)
- [Pipeline Improvements](#pipeline-improvements)

# Overview

<img src="scATAC_diagram.png" width="500">

# Introduction

The scATAC Pipeline was developed by the Broad DSP Pipelines team to process single nucleus ATAC-seq datasets. The pipeline is based on the [SnapATAC pipeline](https://github.com/r3fang/SnapATAC) described by [Fang et al. (2019)](https://www.biorxiv.org/content/10.1101/615179v2.full). Overall, the pipeline uses the python module [SnapTools](https://github.com/r3fang/SnapTools) to align and process paired reads in the form of FASTQ files. It produces an hdf5-structured Snap file that includes a cell-by-bin count matrix at 10 kb resolution. In addition to the Snap file, the final outputs include a GA4GH compliant aligned BAM and QC metrics.

## Quick Start Table

| Pipeline Features | Description | Source |
| ---  |--- | --- |
| Assay Type | Single nucleus ATAC-seq | [Preprint here ](https://www.biorxiv.org/content/biorxiv/early/2019/05/13/615179.full.pdf)
| Overall Workflow  | Generates Snap file with cell x bin matrix at 10 kb resolution | Code available from [GitHub](scATAC.wdl) |
| Workflow Language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Aligner  | BWA | [Li H. and Durbin R., 2009](https://pubmed.ncbi.nlm.nih.gov/19451168/) |                     
| Data Input File Format | File format in which sequencing data is provided | Paired-end FASTQs with cell barcodes appended to read names (read barcode demultiplexing section [here](https://github.com/r3fang/SnapATAC/wiki/FAQs#whatissnap)) |                     
| Data Output File Format | File formats in which scATAC output is provided | [BAM](http://samtools.github.io/hts-specs/), [Snap](https://github.com/r3fang/SnapATAC/wiki/FAQs#whatissnap) |

# Set-up
## Workflow Installation and Requirements
The [scATAC workflow](scATAC.wdl) is written in the Workflow Description Language WDL and can be downloaded by cloning the GitHub [Skylab repository](https://github.com/HumanCellAtlas/skylab). The workflow can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH compliant, flexible workflow management system that supports multiple computing platforms. For the latest workflow version and release notes, please see the scATAC [changelog](scATAC.changelog.md). 

## Pipeline Inputs
The pipeline inputs are detailed in the table below. You can test the workflow by using the [https://github.com/broadinstitute/warp/tree/master/pipelines/skylab/snap-atac/example_inputs/human_example.json](human_example.json) example configuration file. 

| Input name | Input type | Description |
| --- | --- | --- |
| input_fastq1 | File | FASTQ file of the first reads (R1) |
| input_fastq2 | File | FASTQ file of the second reads (R2) |
| input_reference | File | Reference bundle that is generated with bwa-mk-index-wdl found [here](https://github.com/HumanCellAtlas/skylab/blob/master/library/accessory_workflows/build_bwa_reference/bwa-mk-index.wdl)|
| genome_name | String | Name of the genomic reference (name that precedes the “.tar” in the input_reference) |
| output_bam  | String  | Name for the output BAM |

## Input File Preparation

### R1 and R2 FASTQ Preparation
The scATAC workflow requires paired reads in the form FASTQ files with the cell barcodes appended to the readnames. A description of the barcode demultiplexing can be found on the SnapATAC documentation (see barcode demultiplexing section [here](https://github.com/r3fang/SnapATAC/wiki/FAQs#CEMBA_snap)). The full cell barcode must form the first part of the read name (for both R1 and R2 files) and be separated from the rest of the line by a colon. You can find an example python code to perform demultiplexing in the [SnapTools documentation here](https://github.com/r3fang/SnapTools/blob/master/snaptools/dex_fastq.py). The codeblock below demonstrates the correct format. 

```
@CAGTTGCACGTATAGAACAAGGATAGGATAAC:7001113:915:HJ535BCX2:1:1106:1139:1926 1:N:0:0
ACCCTCCGTGTGCCAGGAGATACCATGAATATGCCATAGAACCTGTCTCT
+
DDDDDIIIIIIIIIIIIIIHHIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

### Input_reference Preparation
The input_reference is a BWA compatible reference bundle in TAR file format. You can create this BWA reference using the accessory workflow  [here](https://github.com/HumanCellAtlas/skylab/blob/master/library/accessory_workflows/build_bwa_reference/bwa-mk-index.wdl).


# Workflow Tasks and Tools

The [scATAC workflow](scATAC.wdl) is divided into multiple tasks which are described in the table below. The table also links to the Docker Image for each task and to the documentation or code for the relevant software tool parameters.

| Task | Task Description | Tool Docker Image | Parameter Descriptions or Code |
|--- | --- | --- | --- |
| AlignPairedEnd | Align the modified FASTQ files to the genome | [snaptools:0.0.1](https://github.com/HumanCellAtlas/skylab/blob/master/docker/snaptools/Dockerfile) | [SnapTools documentation](https://github.com/r3fang/SnapTools) |
| SnapPre | Initial generation of snap file | [snaptools:0.0.1](https://github.com/HumanCellAtlas/skylab/blob/master/docker/snaptools/Dockerfile) | [SnapTools documentation](https://github.com/r3fang/SnapTools) |
| SnapCellByBin | Binning of data by genomic bins | [snaptools:0.0.1](https://github.com/HumanCellAtlas/skylab/blob/master/docker/snaptools/Dockerfile) | [SnapTools documentation](https://github.com/r3fang/SnapTools) |
| MakeCompliantBAM | Generation of a GA4GH compliant BAM | [snaptools:0.0.1](https://github.com/HumanCellAtlas/skylab/blob/master/docker/snaptools/Dockerfile) | [Code](https://github.com/HumanCellAtlas/skylab/blob/master/docker/snaptools/makeCompliantBAM.py) |
| BreakoutSnap | Extraction of tables from snap file into text format (for testing and user availability) | [snap-breakout:0.0.1](https://github.com/HumanCellAtlas/skylab/tree/master/docker/snap-breakout) | [Code](https://github.com/HumanCellAtlas/skylab/blob/master/docker/snap-breakout/breakoutSnap.py) |

## Task Summary 

### AlignPairedEnd
The AlignPairedEnd task takes the barcode demultiplexed FASTQ files and aligns reads to the genome using the BWA aligner. It uses the SnapTools min_cov parameter to set the minimum number of barcodes a fragment requires to be included in the final output. This parameter is set to 0. The final task output is an aligned BAM file. 

### SnapPre

The SnapPre task uses SnapTools to perform preprocessing and filtering on the aligned BAM. The task outputs are a Snap file and QC metrics. The table below details the filtering parameters for the task.

#### Filtering Parameters
| Parameter | Description | Value |
| --- | --- | --- |
| --min-mapq | Fragments with mappability less than value will be filtered | 30 |
| --min-flen | Fragments of length shorter than min_flen will be filtered | 0 |
| --max-flen | Fragments of length bigger than min_flen will be filtered | 1000 |
| --keep-chrm | Boolean variable indicates whether to keep reads mapped to chrM | TRUE |
| --keep-single | Boolean variable indicates whether to keep single-end reads | TRUE |
| --keep-secondary | Boolean variable indicates whether to keep secondary alignments | FALSE
| --max-num | Max number of barcodes to be stored. Based on the coverage, top max_barcode barcodes are selected and stored | 1000000 |
| --min-cov | Fragments with less than min-cov number of barcodes will be filtered | 100 |


### SnapCellByBin

The SnapCellByBin task uses the Snap file to create cell-by-bin count matrices in which a “1” represents a bin with an accessible region of the genome and a “0” represents an inaccessible region. The bin_size_list is set to 10,000 bp. 

### MakeCompliantBAM

The MakeCompliantBAM task uses a [custom python script (here)](https://github.com/HumanCellAtlas/skylab/tree/master/docker/snaptools) to make a GA4GH compliant BAM by moving the cellular barcodes in the read names to the CB tag. 

### BreakoutSnap
The BreakoutSnap task extracts data from the Snap file and exports it to individual CSVs. These CSV outputs are listed in the table in the Outputs section below. The code is available [here](https://github.com/HumanCellAtlas/skylab/blob/master/docker/snap-breakout/breakoutSnap.py). 

# Outputs 

The main outputs of the scATAC workflow are the Snap file, Snap QC metrics, and the GA4GH compliant BAM file. All files with the prefix “breakout” are CSV files containing individual pieces of data from the Snap. The sessions for the Snap file are described in the [SnapTools documentation](https://github.com/r3fang/SnapTools). Additionally, you can read detailed information on the [Snap file fields for each session](https://github.com/r3fang/SnapTools/blob/master/docs/snap_format.docx) (select "View Raw").

| Output File Name | Description |
| --- | --- |
| output_snap_qc | Quality control file corresponding to the snap file |
| output_snap | Output snap file (in hdf5 container format) |
| output_aligned_bam  | Output BAM file, compliant with GA4GH standards |
| breakout_barcodes | Text file containing the FM ('Fragment session') barcodeLen and barcodePos fields  |
| breakout_fragments | Text file containing the FM ('Fragments session') fragChrom, fragLen, and fragStart fields |
| breakout_binCoordinates | Text file with the AM session ('Cell x bin accessibility' matrix) binChrom and binStart fields |
| breakout_binCounts  | Text file with the AM session ('Cell x bin accessibility' matrix) idx, idy, and count fields |
| breakout_barcodesSection  | Text file with the data from the BD session ('Barcode session' table) |

#### Snap QC Metrics
The following table details the metrics available in the output_snap_qc file.

| QC Metric | Abbreviation |
| --- | --- |
| Total number of unique barcodes | No abbreviation |
| Total number of fragments | TN |
| Total number of uniquely mapped | UM | 
| Total number of single ends | SE |
| Total number of secondary alignments | SA |
| Total number of paired ends | PE |
| Total number of proper paired | PP |
| Total number of proper frag len | PL |
| Total number of usable fragments | US |
| Total number of unique fragments | UQ |
| Total number of chrM fragments | CM |

# Running on Terra

[Terra](https://app.terra.bio/) is a public, cloud-based platform for biomedical research. If you would like to try the scATAC workflow (previously named "snap-atac") in Terra, you can import the most recent version from the [Broad Methods Repository](https://portal.firecloud.org/?return=terra#methods/snap-atac-v1_0/snap-atac-v1_0/2) (Google login required). Additionally, there is a public [SnapATAC_Pipeline workspace](https://app.terra.bio/#workspaces/brain-initiative-bcdc/SnapATAC_Pipeline) preloaded with the scATAC workflow and downsampled data.  

# Versioning
All scATAC workflow releases are documented in the [scATAC changelog](scATAC.changelog.md).

# Pipeline Improvements
Please help us make our tools better by contacting [Kylee Degatano](mailto:kdegatano@broadinstitute.org) for pipeline-related suggestions or questions.
