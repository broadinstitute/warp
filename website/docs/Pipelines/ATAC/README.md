---
sidebar_position: 1
slug: /Pipelines/ATAC/README

---

# ATAC Overview

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [2.5.1](https://github.com/broadinstitute/warp/releases) | November, 2024 | WARP Pipelines | Please [file an issue in WARP](https://github.com/broadinstitute/warp/issues). |


## Introduction to the ATAC workflow
ATAC is an open-source, cloud-optimized pipeline developed in collaboration with members of the [BRAIN Initiative](https://braininitiative.nih.gov/) (BICCN and [BICAN](https://brainblog.nih.gov/brain-blog/brain-issues-suite-funding-opportunities-advance-brain-cell-atlases-through-centers) Sequencing Working Group) and [SCORCH](https://nida.nih.gov/about-nida/organization/divisions/division-neuroscience-behavior-dnb/basic-research-hiv-substance-use-disorder/scorch-program) (see [Acknowledgements](#acknowledgements) below). It supports the processing of 10x single-nucleus data generated with 10x Multiome [ATAC-seq (Assay for Transposase-Accessible Chromatin)](https://www.10xgenomics.com/products/single-cell-multiome-atac-plus-gene-expression), a technique used in molecular biology to assess genome-wide chromatin accessibility. 

This workflow is the ATAC component of the [Mutiome wrapper workflow](../Multiome_Pipeline/README). It corrects cell barcodes (CBs), aligns reads to the genome, and produces a fragment file as well as [per barcode metrics](../ATAC/count-matrix-overview.md) and [library-level metrics](../ATAC/library-metrics.md).


## Quickstart table
The following table provides a quick glance at the ATAC pipeline features:

| Pipeline features | Description | Source |
|--- | --- | --- |
| Assay type | 10x single cell or single nucleus ATAC | [10x Genomics](https://www.10xgenomics.com)
| Overall workflow  | Barcode correction, read alignment, and fragment quantification | Code available from [GitHub](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/atac/atac.wdl)
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Genomic Reference Sequence | GRCh38 human genome primary sequence | GENCODE |
| Aligner | bwa-mem2 | [Li H. and Durbin R., 2009](http://www.ncbi.nlm.nih.gov/pubmed/19451168) |
| Fragment quantification | SnapATAC2 | [Zhang, K. et al., 2021](https://pubmed.ncbi.nlm.nih.gov/34774128/)
| Data input file format | File format in which sequencing data is provided | [FASTQ](https://academic.oup.com/nar/article/38/6/1767/3112533) |
| Data output file format | File formats in which ATAC output is provided | TSV, h5ad, BAM |
| Library-level metrics | The [ATAC](../ATAC/README.md) pipeline uses [SnapATAC2](https://github.com/kaizhang/SnapATAC2) to generate library-level metrics in CSV format. | [Library-level metrics](../ATAC/library-metrics.md) |

## Set-up

### ATAC installation

To download the latest ATAC release, see the release tags prefixed with "Multiome" on the WARP [releases page](https://github.com/broadinstitute/warp/releases). All ATAC pipeline releases are documented in the [ATAC changelog](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/multiome/atac.changelog.md). 

To discover and search releases, use the WARP command-line tool [Wreleaser](https://github.com/broadinstitute/warp/tree/master/wreleaser).

ATAC can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant, flexible workflow management system that supports multiple computing platforms. The workflow can also be run in [Terra](https://app.terra.bio), a cloud-based analysis platform. 

## Input Variables
The following describes the inputs of the ATAC workflow. For more details on how default inputs are set for the Multiome workflow, see the [Multiome overview](../Multiome_Pipeline/README).

| Variable name | Description |
| --- |--- |
| read1_fastq_gzipped | Fastq inputs (array of compressed read 1 FASTQ files). |
| read2_fastq_gzipped | Fastq inputs (array of compressed read 2 FASTQ files containing cellular barcodes). |
| read3_fastq_gzipped | Fastq inputs (array of compressed read 3 FASTQ files). |
| input_id | Output prefix/base name for all intermediate files and pipeline outputs.                                        |
| cloud_provider | String describing the cloud provider that should be used to run the workflow; value should be "gcp" or "azure". | String |
| preindex | Boolean used for paired-tag data and not applicable to ATAC data types; default is set to false. | 
| atac_expected_cells | Number of cells loaded to create the ATAC library; default is set to 3000. |
| tar_bwa_reference | BWA reference (tar file containing reference fasta and corresponding files). |
| num_threads_bwa | Optional integer defining the number of CPUs per node for the BWA-mem alignment task (default: 128). |
| mem_size_bwa | Optional integer defining the memory size for the BWA-mem alignment task in GB (default: 512). |
| cpu_platform_bwa | Optional string defining the CPU platform for the BWA-mem alignment task (default: "Intel Ice Lake"). |
| annotations_gtf | CreateFragmentFile input variable: GTF file for SnapATAC2 to calculate TSS sites of fragment file. |
| chrom_sizes | CreateFragmentFile input variable: Text file containing chrom_sizes for genome build (i.e., hg38) |
| whitelist | Whitelist file for ATAC cellular barcodes. |
| adapter_seq_read1 | TrimAdapters input: Sequence adapter for read 1 fastq. |
| adapter_seq_read3 | TrimAdapters input: Sequence adapter for read 3 fastq. |
| vm_size | String defining the Azure virtual machine family for the workflow (default: "Standard_M128s"). 
| atac_nhash_id | String that represents an optional library aliquot identifier. When used, it is echoed in the h5ad unstructured data. |
| peak_calling | Optional boolean used to determine if the ATAC pipeline should run Peak Calling; default is "false".  | Boolean |

## ATAC tasks and tools

Overall, the ATAC workflow:
1. Identifies optimal parameters for performing CB correction and alignment. 
1. Corrects CBs and partitions FASTQs by CB.
1. Aligns reads.
1. Generates a fragment file.
1. Calculates per cell barcode fragment metrics.

The tools each ATAC task employs are detailed in the table below. 

To see specific tool parameters, select the task WDL link in the table; then view the `command {}` section of the task in the WDL script. To view or use the exact tool software, see the task's Docker image which is specified in the task WDL `# runtime values` section as `String docker =`.

| Task name and WDL link | Tool | Software | Description | 
| --- | --- | --- | ------------------------------------ | 
| [GetNumSplits](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/atac/atac.wdl) | Bash | Bash | Uses the virtual machine type to determine the optimal number of FASTQ files for performing the BWA-mem alignment step. This allows BWA-mem to run in parallel on multiple FASTQ files in the subsequent workflow steps. |
| [FastqProcessing as SplitFastq](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/FastqProcessing.wdl) | fastqprocess | custom | Dynamically selects the correct barcode orientation, corrects cell barcodes, and splits FASTQ files by the optimal number determined in the [GetNumSplits](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/atac/atac.wdl) task. The smaller FASTQ files are grouped by cell barcode with each read having the corrected (CB) and raw barcode (CR) in the read name. |
| [TrimAdapters](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/atac/atac.wdl) | Cutadapt v4.4 | cutadapt | Trims adaptor sequences. |
| [BWAPairedEndAlignment](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/atac/atac.wdl) | bwa-mem2 | mem | Aligns reads from each set of partitioned FASTQ files to the genome and outputs a BAM with ATAC barcodes in the CB:Z tag. |
| [CreateFragmentFile](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/atac/atac.wdl) | make_fragment_file, import_data | SnapATAC2 | Generates a fragment file from the final aligned BAM and outputs per barcode quality metrics in h5ad. A detailed list of these metrics is found in the [ATAC Count Matrix Overview](./count-matrix-overview.md). This task is nondeterministic.|
| [PeakCalling](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/atac/atac.wdl) | macs3 | SnapATAC2 | Generates two h5ad files (`cellbybin.h5ad` and `cellbypeak.h5ad`) from the CreateFragmentFile h5ad output file (`metrics.h5ad`). The `cellbybin.h5ad` contains the peak called per cluster in the macs3 unstructured metadata and `cellbypeak.h5ad` contains the merged peaks and the count matrix of peaks per fragment. A detailed list of these metrics is found in the [ATAC Count Matrix Overview](./count-matrix-overview.md). This task is nondeterministic.|

## Output variables

| Output variable name | Filename, if applicable | Output format and description |
|--- | --- | --- | 
| bam_aligned_output | `<input_id>`.bam | BAM containing aligned reads from ATAC workflow. |
| fragment_file | `<input_id>`.fragments.sorted.tsv.gz | Bgzipped TSV containing fragment start and stop coordinates per barcode. In order, the columns are "Chromosome", "Start", "Stop", "ATAC Barcode", and "Number Reads". | 
| snap_metrics | `<input_id`.metrics.h5ad | h5ad (Anndata) containing per barcode metrics from [SnapATAC2](https://github.com/kaizhang/SnapATAC2). A detailed list of these metrics is found in the [ATAC Count Matrix Overview](./count-matrix-overview.md). |
 library_metrics | `<input_id>`_`<atac_nhash_id>`_library_metrics.csv | CSV file containing library-level metrics. Read more in the [Library Metrics Overview](library-metrics.md) |
| snap_cellbybin | `<input_id>`.cellbybin.h5ad | h5ad (Anndata) containing peaks (called by MACS3) per cluster. [SnapATAC2](https://github.com/kaizhang/SnapATAC2). |
| snap_cellbypeak | `<input_id>`.cellbypeak.h5ad | h5ad (Anndata) containing merged peaks (called by MACS3) per cluster and count matrix of insertion sites per peak and cell. [SnapATAC2](https://github.com/kaizhang/SnapATAC2).|

## Versioning and testing

All ATAC pipeline releases are documented in the [ATAC changelog](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/multiome/atac.changelog.md) and tested using [plumbing and scientific test data](https://github.com/broadinstitute/warp/tree/master/pipelines/skylab/multiome/test_inputs). To learn more about WARP pipeline testing, see [Testing Pipelines](https://broadinstitute.github.io/warp/docs/About_WARP/TestingPipelines).

## Citing the ATAC Pipeline

If you use the ATAC Pipeline in your research, please identify the pipeline in your methods section using the [ATAC SciCrunch resource identifier](https://scicrunch.org/resources/data/record/nlx_144509-1/SCR_024656/resolver?q=SCR_024656%2A&l=SCR_024656%2A&i=rrid:scr_024656).

* Ex: *ATAC Pipeline (RRID:SCR_024656)*

When citing WARP, please use the following:

Degatano, Kylee, Aseel Awdeh, Robert Sidney Cox III, Wes Dingman, George Grant, Farzaneh Khajouei, Elizabeth Kiernan, et al. 2025. "Warp Analysis Research Pipelines: Cloud-Optimized Workflows for Biological Data Processing and Reproducible Analysis." _Bioinformatics (Oxford, England)_, September, https://doi.org/10.1093/bioinformatics/btaf494 .


## Acknowledgements

We are immensely grateful to the members of the BRAIN Initiative (BICAN Sequencing Working Group) and SCORCH for their invaluable and exceptional contributions to this pipeline. Our heartfelt appreciation goes to Alex Dobin, Aparna Bhaduri, Alec Wysoker, Anish Chakka, Brian Herb, Daofeng Li, Fenna Krienen, Guo-Long Zuo, Jeff Goldy, Kai Zhang, Khalid Shakir, Bo Li, Mariano Gabitto, Michael DeBerardine, Mengyi Song, Melissa Goldman, Nelson Johansen, James Nemesh, and Theresa Hodges for their unwavering dedication and remarkable efforts.

