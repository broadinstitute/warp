| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [optimus_v4.0.2](https://github.com/broadinstitute/warp/releases) | September 14, 2020 | [Elizabeth Kiernan](mailto:ekiernan@broadinstitute.org) | Please file GitHub issues in warp or contact [Kylee Degatano](mailto:kdegatano@broadinstitute.org) |


# Table of Contents
- [Optimus Pipeline Overview](#optimus-pipeline-overview)
  * [Introduction to the Optimus Workflow](#introduction-to-the-optimus-workflow)
  * [Quick Start Table](#quick-start-table)
- [Set-up](#set-up)
  * [Optimus Installation and Requirements](#optimus-installation-and-requirements)
  * [Inputs](#inputs)
    + [Sample Data Input](#sample-data-input)
    + [Additional Reference Inputs](#additional-reference-inputs)
    + [Sample Inputs for Analyses in a Terra Workspace](#sample-inputs-for-analyses-in-a-terra-workspace)
- [Optimus Tasks and Tools](#optimus-tasks-and-tools)
  * [Optimus Task Summary](#optimus-task-summary)
    + [1. Converting R2 FASTQ File to UBAM](#1-converting-r2-fastq-file-to-ubam)
    + [2. Correcting and Attaching Cell Barcodes](#2-correcting-and-attaching-cell-barcodes)
    + [3. Alignment](#3-alignment)
    + [4. Gene Annotation](#4-gene-annotation)
    + [5. UMI Correction](#5-umi-correction)
    + [6. Summary Metric Calculation](#6-summary-metric-calculation)
    + [7. Expression Matrix Construction](#7-expression-matrix-construction)
    + [8. Identification of Empty Droplets](#8-identification-of-empty-droplets)
    + [9. Outputs](#9-outputs)
- [Versioning](#versioning)
- [Pipeline Improvements](#have-suggestions)
- [FAQs](#faqs)


# Optimus Pipeline Overview
![Diagram](documentation/Optimus_diagram.png)

## Introduction to the Optimus Workflow

Optimus is a pipeline developed by the Data Coordination Platform (DCP) of the [Human Cell Atlas (HCA) Project](https://data.humancellatlas.org/) that supports processing of any 3' single-cell and single-nuclei expression data generated with the [10x Genomic v2 or v3 assay](https://www.10xgenomics.com/solutions/single-cell/). It is an alignment and transcriptome quantification pipeline that corrects cell barcodes, aligns reads to the genome, corrects Unique Molecular Identifiers (UMIs), generates an expression matrix in a UMI-aware manner, calculates summary metrics for genes and cells, detects empty droplets, returns read outputs in BAM format, and returns cell gene expression in numpy matrix and Loom file formats. Special care is taken to keep all reads that may be useful to the downstream user, such as unaligned reads or reads with uncorrectable barcodes. This design provides flexibility to the downstream user and allows for alternative filtering or leveraging the data for novel methodological development.

Optimus has been validated for analyzing both [human](https://github.com/broadinstitute/warp/blob/master/benchmarking/v1_Apr2019/optimus_report.rst) and [mouse](https://docs.google.com/document/d/1_3oO0ZQSrwEoe6D3GgKdSmAQ9qkzH_7wrE7x6_deL10/edit) single-cell or single-nuclei data sets. See the [human single-cell validation](https://docs.google.com/document/d/158ba_xQM9AYyu8VcLWsIvSoEYps6PQhgddTr9H0BFmY/edit) or the [single-nuclei](https://docs.google.com/document/d/1rv2M7vfpOzIOsMnMfNyKB4HV18lQ9dnOGHK2tPikiH0/edit) validation reports. 

| Want to use the Optimus Pipeline for your publication? |
|---|
| *Check out the [Optimus Publication Methods](optimus.methods.md) to get started!* |

## Quick Start Table

| Pipeline Features | Description | Source |
|-------------------|---------------------------------------------------------------|-----------------------|
| Assay Type | 10x Single Cell or Single Nuclei Expression (v2 and v3) | [10x Genomics](https://www.10xgenomics.com)
| Overall Workflow  | Quality control module and transcriptome quantification module | Code available from [GitHub](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/Optimus.wdl) |
| Workflow Language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Genomic Reference Sequence| GRCh38 human genome primary sequence and M21 (GRCm38.p6) mouse genome primary sequence | GENCODE [Human](https://www.gencodegenes.org/human/release_27.html) and [Mouse](https://www.gencodegenes.org/mouse/release_M21.html) 
| Transcriptomic Reference Annotation | V27 GENCODE human transcriptome and M21 mouse transcriptome | GENCODE [Human](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz) and [Mouse](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.gff3.gz) |
| Aligner  | STAR (v.2.5.3a) | [Dobin, et al.,2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) |
| Transcript Quantification | Utilities for processing large-scale single cell datasets | [sctools](https://github.com/HumanCellAtlas/sctools)                          
| Data Input File Format | File format in which sequencing data is provided | [FASTQ](https://academic.oup.com/nar/article/38/6/1767/3112533) |                     
| Data Output File Format | File formats in which Optimus output is provided | [BAM](http://samtools.github.io/hts-specs/), Python numpy arrays (internal), Loom (generated with [Loompy v.3.0.6)](http://loompy.org/) |

# Set-up

## Optimus Installation and Requirements
The Optimus pipeline code can be downloaded by cloning the GitHub repository [warp](https://github.com/broadinstitute/warp/). For the latest release of Optimus, please see the release tags prefixed with "Optimus" [here](https://github.com/broadinstitute/warp/releases). 

Optimus can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH compliant, flexible workflow management system that supports multiple computing platforms. Optimus can also be run in [Terra](https://app.terra.bio), a cloud-based analysis platform. The Terra [Optimus Featured Workspace](https://app.terra.bio/#workspaces/featured-workspaces-hca/HCA_Optimus_Pipeline) contains the Optimus workflow, workflow configurations, required reference data and other inputs, and example testing data.

## Inputs

Optimus pipeline inputs are detailed in JSON format configuration files. There are four example configuration files available if you are interested in running the pipeline: 
*  [human_v2_example](./example_inputs/human_v2_example.json): An example human 10x v2 single-cell dataset
*  [human_v3_example](./example_inputs/human_v3_example.json): An example human 10x v3 single-cell dataset
*  [mouse_v2_example](./example_inputs/mouse_v2_example.json): An example mouse 10x v2 single-cell dataset
*  [mouse_v2_snRNA_example](./example_inputs/mouse_v2_snRNA_example.json): An example mouse v2 single-nuclei dataset

Additionally, there are multiple sample datasets available in the [test_optimus_full_datasets](./example_inputs/test_optimus_full_datasets) folder. Please note that unlike the example configuration files above, the configuration files in this folder may not reflect updated Optimus parameters. However, you can still access the FASTQ files for each dataset at the Google bucket locations listed in the dataset configuration files. 

### Sample Data Input

Each 10x v2 and v3 3’ sequencing experiment generates triplets of FASTQ files for any given sample:  

1. Forward reads (r1_fastq) containing the unique molecular identifier (UMI) and cell barcode sequences
2. Reverse reads (r2_fastq) containing the alignable genomic information from the mRNA transcript 
3. Index FASTQ (i1_fastq) containing the sample barcodes, when provided by the sequencing facility

Note: Optimus is currently a single sample pipeline, but can take in multiple sets of FASTQs for a sample that has been split over multiple lanes of sequencing. For an example configuration file with multiple lanes, please see the [mouse_v2_example.json](mouse_v2_example.json). Additionally, Optimus does not support demultiplexing even though it accepts index FASTQ files. 

### Additional Reference Inputs

The JSON file also contains metadata for the reference information in the following table:

| Parameter Name | Description | Optional Strings (when applicable) |
| --- | --- | --- |
| Whitelist | Cloud path to list of known cell barcodes from [10x genomics](https://www.10xgenomics.com/) that corresponds to the v2 or v3 chemistry | NA |
| Tar_star_reference | Cloud path to TAR file containing a species-specific reference genome and gtf; it is generated using the [BuildIndices.wdl](../build_indices/BuildIndices.wdl) | NA |
| input_id | Unique name describing the biological sample or replicate that corresponds with the original FASTQ files. This can be any string, but if possible, we recommend it matches the sample metadata | NA | 
| input_name | Optional string that can be used to further identify the original biological sample | NA |
| Annotations_gtf | Cloud path to GTF containing gene annotations used for gene tagging (must match GTF in STAR reference) | NA | 
| Chemistry | Optional string description of whether data was generated with 10x v2 or v3 chemistry. Optimus validates this string. If the string does not match one of the optional strings, the pipeline will fail. You can remove the checks by setting "force_no_check = true" in the input JSON | "tenX_v2" (default) or "tenX_v3" |
| Counting_mode | String description of whether data is single-cell or single-nuclei | "sc_rna" or "sn_rna" |
| Output_bam_basename | Optional string used for the output BAM file basename; the default is input_id | NA |


### Sample Inputs for Analyses in a Terra Workspace

The Optimus pipeline is currently available on the cloud-based platform Terra. If you have a Terra account, you can access the Featured Workspace using this address: [https://app.terra.bio/#workspaces/featured-workspaces-hca/HCA_Optimus_Pipeline](https://app.terra.bio/#workspaces/featured-workspaces-hca/HCA_Optimus_Pipeline). The workspace is preloaded with instructions and sample data. For more information on using the Terra platform, please view the [Support Center](https://support.terra.bio/hc/en-us).  

# Optimus Tasks and Tools

The [Optimus.wdl](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/Optimus.wdl) in the pipelines/optimus folder of the WARP repository implements the workflow by importing individual modules ("tasks" written in  WDL script) from the WARP [tasks](https://github.com/broadinstitute/warp/blob/master/tasks/skylab) folder.

## Optimus Task Summary

Here we describe the tasks of the Optimus pipeline; [the code](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/Optimus.wdl) and [library of tasks](https://github.com/broadinstitute/warp/blob/master/tasks/skylab) are available through GitHub.

Overall, the workflow:
1. Converts R2 FASTQ file (containing alignable genomic information) to an unaligned BAM (UBAM)
2. Corrects and attaches 10x Barcodes using the R1 FASTQ file 
3. Aligns reads to the genome with STAR v.2.5.3a
4. Annotates genes with aligned reads
5. Corrects UMIs
6. Calculates summary metrics
7. Produces a UMI-aware expression matrix
8. Detects empty droplets
9. Returns a GA4GH compliant BAM and an expression matrix in Loom formats

The tools each Optimus task employs are detailed in the table below. If you are looking for the parameters for each task/tool, please click on the task link and see the `command {}` section of the task WDL script. The task's Docker image is specified in the task WDL `# runtime values` section as ``` String docker = ```. 

| Task | Tool | 
| --- | --- |
| [FastqProcessing](/tasks/skylab/FastqProcessing.wdl) | [sctools](https://sctools.readthedocs.io/en/latest/sctools.html) | 
| [StarAlignBamSingleEnd](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/StarAlignBamSingleEnd.wdl) |	[STAR](https://github.com/alexdobin/STAR) |
| [TagGeneExon](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/TagGeneExon.wdl) |	[Drop-seq](https://github.com/broadinstitute/Drop-seq) |
| [UmiCorrection](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/UmiCorrection.wdl) |	[Umi-tools](https://github.com/CGATOxford/UMI-tools) |
| [SequenceDataWithMoleculeTagMetrics](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/SequenceDataWithMoleculeTagMetrics.wdl) |	[sctools](https://sctools.readthedocs.io/en/latest/sctools.html) |
| [RunEmptyDrops](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/RunEmptyDrops.wdl) |	[dropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html) |
| [CreateCountMatrix](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/CreateCountMatrix.wdl) |	[Drop-seq](https://github.com/broadinstitute/Drop-seq) and [sctools](https://sctools.readthedocs.io/en/latest/sctools.html)|
| [FastqToUBAM](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/FastqToUBam.wdl)	| [picard](https://github.com/broadinstitute/picard) |
| [TagSortBam](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/TagSortBam.wdl) |	[sctools](https://sctools.readthedocs.io/en/latest/sctools.html) |
| [Picard](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/Picard.wdl)	| [picard](https://github.com/broadinstitute/picard) |

The Optimus pipeline takes special care to flag but avoid the removal of reads that are not aligned or that do not contain recognizable barcodes. This design (which differs from many pipelines currently available) allows the use of the entire dataset by those who may want to use alternative filtering or leverage the data for methodological development associated with the data processing. More information about the different tags used to flag the data can be found in the [Bam_tags documentation](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/documentation/Bam_tags.md).

### 1. Converting R2 FASTQ File to UBAM

Unlike FASTQ files, BAM files enable researchers to keep track of important metadata throughout all data processing steps. The first step of Optimus is to [convert](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/FastqProcessing.wdl) the R2 FASTQ file, containing the alignable genomic information, to an unaligned BAM (UBAM) file.

### 2. Correcting and Attaching Cell Barcodes

Although the function of the cell barcodes is to identify unique cells, barcode errors can arise during sequencing (such as incorporation of the barcode into contaminating DNA or sequencing and PCR errors), making it difficult to distinguish unique cells from artifactual appearances of the barcode. The [FastqProcessing](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/FastqProcessing.wdl) task uses [sctools](https://github.com/HumanCellAtlas/sctools) to evaluate barcode errors by comparing the R1 FASTQ sequences against a whitelist of known barcode sequences. The task then appends the UMI and cell barcode sequences from the R1 FASTQ to the UBAM sequence as tags [(see the Bam_tags documentation for details](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/documentation/Bam_tags.md)). 

The output is a UBAM file containing the reads with corrected barcodes, including barcodes that came within one edit distance ([Levenshtein distance](http://www.levenshtein.net/)) of matching the whitelist of barcode sequences and were corrected by this tool. Correct barcodes are assigned a “CB” tag. Uncorrectable barcodes (with more than one error) are preserved and given a “CR” (Cell barcode Raw) tag. Cell barcode quality scores are also preserved in the file under the “CY” tag.


### 3. Alignment

Optimus uses the [STAR alignment](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/StarAlignBamSingleEnd.wdl) task to map barcoded reads in the UBAM file to the genome primary assembly reference (see table above for version information). This task uses STAR (Spliced Transcripts Alignment to a Reference) a standard, splice-aware, RNA-seq alignment tool [(Dobin, et al., 2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/). The default soft-clipping is turned on for this alignment. 


### 4. Gene Annotation

The [TagGeneExon](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/TagGeneExon.wdl) task uses [Drop-seq tools](https://github.com/broadinstitute/Drop-seq) to [annotate each read](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/documentation/Bam_tags.md) with the type of sequence to which it aligns. These annotations vary depending on the counting_mode ("sc_rna" or "sn_rna") specified in the workflow.

**Single-cell RNA-seq:**

The TagGeneExon task calls Drop-seq tools v1.12 to make annotations. These annotations include INTERGENIC, INTRONIC, UTR and CODING (EXONIC), and are stored using the 'XF' BAM tag. In cases where the gene corresponds to an exon or UTR, the name of the gene that overlaps the alignment is associated with the read and stored using the GE BAM tag. 

**Single-nuclei RNA-seq:**

The TagGeneExon task calls Drop-seq tools v2.3.0 to make annotations. These annotations include INTERGENIC, INTRONIC, UTR and CODING (EXONIC), and are stored using the 'XF' BAM tag (see the [Bam_tags documentation](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/documentation/Bam_tags.md)). In cases where the gene corresponds to an exon, UTR, or intron, the name of the gene that overlaps the alignment is associated with the read and stored using the 'GE' BAM tag. 

All tags are detailed in the pipeline's [BAM_tag documentation](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/documentation/Bam_tags.md).

### 5. UMI Correction

UMIs are designed to distinguish unique transcripts present in the cell at lysis from those arising from PCR amplification of these same transcripts. But, like cell barcodes, UMIs can also be incorrectly sequenced or amplified. The [UmiCorrection](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/UmiCorrection.wdl) task uses [Umi-tools v.0.0.1](https://pypi.org/project/umi-tools/0.0.1/) to apply a network-based, "directional" correction method ([Smith, et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5340976/)) to account for such errors. This task makes UMI corrections to alignments made with the 'GE' tag. This step will add a 'UB' tag for UMI-corrected barcodes. 

### 6. Summary Metric Calculation

The [Metrics](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/SequenceDataWithMoleculeTagMetrics.wdl) task uses [sctools](https://github.com/HumanCellAtlas/sctools) to calculate summary metrics which help assess the quality of the data output each time this pipeline is run. These metrics are included in the Loom output file. A detailed list of these metrics is found in the [Loom_schema documentation](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/documentation/Loom_schema.md).

### 7. Expression Matrix Construction

The Optimus [CreateCountMatrix](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/CreateCountMatrix.wdl) task (imported as "Count") evaluates every read in the BAM file and creates a UMI-aware expression matrix using [sctools](https://github.com/HumanCellAtlas/sctools). This matrix contains the number of molecules that were observed for each cell barcode and for each gene. The task discards any read that maps to more than one gene, and counts any remaining reads provided the triplet of cell barcode, molecule barcode, and gene name is unique, indicating the read originates from a single transcript present at the time of cell lysis. To correctly specific the gene name tag, this task will look for the 'GE' tag.

### 8. Identification of Empty Droplets

Empty droplets are lipid droplets that did not encapsulate a cell during 10x sequencing, but instead acquired cell-free RNA (secreted RNA or RNA released during cell lysis) from the solution in which the cells resided ([Lun, et al., 2018](https://www.ncbi.nlm.nih.gov/pubmed/?term=30902100). This ambient RNA can serve as a substrate for reverse transcription, leading to a small number of background reads. The Optimus pipeline calls the [RunEmptyDrops](https://github.com/broadinstitute/warp/blob/master/tasks/skylab/RunEmptyDrops.wdl) task which uses the [dropletUtils v.0.1.1](http://bioconductor.org/packages/release/bioc/html/DropletUtils.html) R package to flag cell barcodes that represent empty droplets rather than cells. A cell will be flagged if it contains fewer than 100 molecules. These metrics are stored in the output Loom file. Details of all the metrics included in the final output files can be found in the [Loom_schema documentation](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/documentation/Loom_schema.md). 

| Warning: RunEmptyDrops output not included for single-nuclei data |
| --- |
| EmptyDrops relies on a visual knee point inflection (described in [Lun et al. (2019)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y)) to differentiate ambient-like cells from empty droplets. If snRNAseq data does not produce a visual knee point inflection, EmptyDrops drops may not accurately distinguish ambient-like cells. When the Optimus workflow counting_mode = sn_rna, it does not include the RunEmptyDrops output in the final Zarr or Loom | 

### 9. Outputs

Output files of the pipeline include:

1. Cell x Gene unnormalized, but UMI-corrected, expression matrices
2. Unfiltered, sorted BAM file with barcode and downstream analysis [tags](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/documentation/Bam_tags.md)
3. Cell metadata, including cell metrics
4. Gene metadata, including gene metrics

The following table lists the output files produced from the pipeline. For samples that have sequenced over multiple lanes, the pipeline will output one merged version of each listed file. 

| Output Name | Filename, if applicable | Output Type |Output Format |
| ------ |------ | ------ | ------ | 
| pipeline_version | | Version of the processing pipeline run on this data | String | 
| bam | <input_id>.bam | Aligned BAM | BAM |
| matrix_row_index | sparse_counts_row_index.npy | Index of cells in expression matrix | Numpy array index |
| matrix_col_index | sparse_counts_col_index.npy | Index of genes in expression matrix | Numpy array index | 
| cell_metrics | merged-cell-metrics.csv.gz | cell metrics | compressed csv | Matrix of metrics by cells | 
| gene_metrics | merged-gene-metrics.csv.gz | gene metrics | compressed csv | Matrix of metrics by genes | 
| loom_output_file | <input_id>.loom | Loom | Loom | Loom file with expression data and metadata | N/A |


The Loom is the default output. See the [create_loom_optimus.py](https://github.com/broadinstitute/warp/blob/master/dockers/skylab/loom-output/create_loom_optimus.py) for the detailed code. The final Loom output contains the unnormalized (unfiltered), UMI-corrected expression matrices, as well as the gene and cell metrics detailed in the [Loom_schema documentation](https://github.com/broadinstitute/warp/blob/master/pipelines/skylab/optimus/documentation/Loom_schema.md). 

| Zarr Array Deprecation Notice June 2020 |
| --- |
| Please note that we have deprecated the previously used Zarr array output. The pipeline now uses the Loom file format as the default output. |

# Versioning

All Optimus pipeline releases are documented in the [Optimus changelog](Optimus.changelog.md). 

# Have Suggestions? 

Coming soon, we will have a GitHub document dedicated to open issues! In the meantime, please help us make our tools better by contacting [Kylee Degatano](mailto:kdegatano@broadinstitute.org) for pipeline-related suggestions or questions.


# FAQs
<details>
<summary>Can I run Optimus in Terra?</summary>
<br>
Yes! We have a Terra workspace that is preconfigured with the latest Optimus workflow and is preloaded with human and mouse sample data. You can access the workspace at https://app.terra.bio/#workspaces/featured-workspaces-hca/HCA_Optimus_Pipeline. You will need a Google account to set up Terra. Please see Terra Support at https://support.terra.bio/hc/en-us for documents on getting started.
</details>

<details>
<summary>Is the output count matrix filtered or normalized?</summary>
<br>
No, we do not filter. We keep as much data as possible so that the researcher can make their own filtering and normalization choices. We do, however, output some information that may be helpful for filtering, like UMI counts per cell and calls on whether or not a cell is empty from EmptyDrops software. For the EmptyDrops call, a cell will be flagged as possibly empty if it contains fewer than 100 molecules.
</details>

<details>
<summary>How does the workflow change when using the single-cell RNA-seq (counting_mode = 'sc_rna') vs. the single-nuclei (counting_mode = 'sn_rna') parameters? </summary>
<br>
Three Optimus tasks are affected by the counting_mode parameter: TagGeneExon, UMICorrection and CreateCountMatrix. The TagGeneExon task uses different versions of Drop-seq tools depending on the counting_mode parameter. The sc_rna parameter uses v1.12 whereas the sn_rna uses v2.3.0. For the sn_rna parameter, a GE tag is added to intronic reads. For the UMICorrection and CreateCountMatrix tasks, the only difference related to the counting_mode parameters is that the sn_RNA parameter will have a GE tag on intronic reads, which the UMICorrection and CreateCountMatrix will recognize. 
Also note that although the RunEmptyDrops task is unaffected by the sn_rna parameter, the final Zarr and Loom outputs will not include EmptyDrops data for the sn_rna mode.
 </details>

<details>
<summary>Where can I find example Optimus datasets and parameters to test the pipeline? </summary>
<br>
There are four example configuration JSON files available for you to test the pipeline- the <a href="./example_inputs/human_v2_example.json">human_v2_example.json </a>, the <a href="./example_inputs/human_v3_example.json">human_v3_example.json </a>, the <a href="./example_inputs/mouse_v2_example.json">mouse_v2_example.json </a>, and the <a href="./example_inputs/mouse_v2_snRNA_example.json">mouse_v2_snRNA_example.json </a>(see the Inputs section). Each of these configuration files can be run in the Optimus Featured Workspace in Terra at https://app.terra.bio/#workspaces/featured-workspaces-hca/HCA_Optimus_Pipeline, but you should note that the workspace comes preloaded with the same data and configurations. We also have multiple example datasets available in the <a href="./example_inputs/test_optimus_full_datasets/">test_optimus_full_datasets </a> folder. These datasets contain configuration files listing the cloud location of the dataset FASTQ files; however, the configuration files may not be updated with all the workflow parameters. For the most up-to-date configuration examples, see the four example files listed above. 
 </details>
 
 <details>
<summary>What outputs are expected if my sample has been sequenced over multiple lanes? </summary>
<br>
The Optimus pipeline is a single sample pipeline, but it can accept multiple FASTQ files if a sample is sequenced across lanes. In this case, the pipeline will merge the results from each lane into single output files. There will only be one merged file for each output type (i.e one Loom, etc.). If you would like to view an example configuration file for a multi-lane dataset, please see the <a href="./example_inputs/mouse_v2_example.json">mouse_v2_example.json </a>.  Additionally, you can view sample outputs in the Optimus featured workspace on Terra: https://app.terra.bio/#workspaces/featured-workspaces-hca/HCA_Optimus_Pipeline. 
 </details>

 <details>
<summary>How do I find which parameters and Docker images were used for the different tasks (i.e. STAR alignment, emptyDrops, etc.) </summary>
<br>
Parameters are listed in each task WDL. For a list of the tasks, see the table in the <a href="README.md/#optimus-task-summary">Task Summary Section </a>. Select the link for the task of interest and then view the parameters in the task WDL "command {}" section. For the task Docker image, see task WDL "# runtime values" section; the Docker is listed as "String docker =  ". If you want to learn more about all the different parameters available for a software tool, please select the relevant link in the table's "Tool" column. 
 </details>




