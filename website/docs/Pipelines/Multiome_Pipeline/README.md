---
sidebar_position: 1
slug: /Pipelines/Multiome_Pipeline/README
---

# Multiome Overview

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [Multiome v2.2.2](https://github.com/broadinstitute/warp/releases) | October, 2023 | Kaylee Mathews | Please file GitHub issues in warp or contact the [WARP Pipeline Development team](mailto:warp-pipelines-help@broadinstitute.org) |

![Multiome_diagram](./multiome_diagram.png)

## Introduction to the Multiome workflow

Multiome is an open-source, cloud-optimized pipeline developed in collaboration with members of the [BRAIN Initiative](https://braininitiative.nih.gov/) (BICCN and [BICAN](https://brainblog.nih.gov/brain-blog/brain-issues-suite-funding-opportunities-advance-brain-cell-atlases-through-centers) Sequencing Working Group) and [SCORCH](https://nida.nih.gov/about-nida/organization/divisions/division-neuroscience-behavior-dnb/basic-research-hiv-substance-use-disorder/scorch-program) (see [Acknowledgements](#acknowledgements) below). It supports the processing of 10x 3' single-cell and single-nucleus gene expression (GEX) and chromatin accessibility (ATAC) data generated with the [10x Genomics Multiome assay](https://www.10xgenomics.com/products/single-cell-multiome-atac-plus-gene-expression).

The workflow is a wrapper WDL script that calls two subworkflows: the [Optimus workflow](../Optimus_Pipeline/README) for single-cell GEX data and the [ATAC workflow](../ATAC/README) for single-cell ATAC data. 

The GEX component corrects cell barcodes (CBs) and Unique Molecular Identifiers (UMIs), aligns reads to the genome, calculates per-barcode and per-gene quality metrics, and produces a raw cell-by-gene count matrix. 

The ATAC component corrects CBs, aligns reads to the genome, calculates per-barcode quality metrics, and produces a fragment file.

The wrapper WDL is available in the WARP repository (see the [code here](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/multiome/Multiome.wdl)).

## Quickstart table
The following table provides a quick glance at the Multiome pipeline features:

| Pipeline features | Description | Source |
|--- | --- | --- |
| Assay type | 10x single cell or single nucleus gene expression (GEX) and ATAC | [10x Genomics](https://www.10xgenomics.com) |
| Overall workflow  | Barcode correction, read alignment, gene and fragment quanitification |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Genomic Reference Sequence | GRCh38 human genome primary sequence | GENCODE [human reference files](https://www.gencodegenes.org/human/release_43.html)|
| Gene annotation reference (GTF) | Reference containing gene annotations | GENCODE [human GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz) |
| Aligners | STARsolo (GEX), BWA-mem2 (ATAC) | [Kaminow et al. 2021](https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1), [Vasimuddin et al. 2019](https://ieeexplore.ieee.org/document/8820962) |
| Transcript and fragment quantification | STARsolo (GEX), SnapATAC2 (ATAC) | [Kaminow et al. 2021](https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1), [SnapATAC2](https://kzhang.org/SnapATAC2/) |
| Data input file format | File format in which sequencing data is provided | [FASTQ](https://academic.oup.com/nar/article/38/6/1767/3112533) |
| Data output file format | File formats in which Multiome output is provided | [BAM](http://samtools.github.io/hts-specs/) and [h5ad](https://anndata.readthedocs.io/en/latest/) |

## Set-up

### Multiome installation

To download the latest Multiome release, see the release tags prefixed with "Multiome" on the WARP [releases page](https://github.com/broadinstitute/warp/releases). All Multiome pipeline releases are documented in the [Multiome changelog](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/multiome/Multiome.changelog.md). 

To discover and search releases, use the WARP command-line tool [Wreleaser](https://github.com/broadinstitute/warp/tree/master/wreleaser).

If you’re running a Multiome workflow version prior to the latest release, the accompanying documentation for that release may be downloaded with the source code on the WARP [releases page](https://github.com/broadinstitute/warp/releases) (see the source code folder).

Multiome can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH compliant, flexible workflow management system that supports multiple computing platforms. The workflow can also be run in [Terra](https://app.terra.bio), a cloud-based analysis platform. 

## Inputs

| Input name | Description | Type |
| --- | --- | --- |
| input_id | Unique identifier describing the biological sample or replicate that corresponds with the FASTQ files; can be a human-readable name or UUID. | String |
| annotations_gtf | GTF file containing gene annotations used for GEX cell metric calculation and ATAC fragment metrics; must match the GTF used to build the STAR aligner. | File |
| gex_r1_fastq | Array of read 1 FASTQ files representing a single GEX 10x library. | Array[File] |
| gex_r2_fastq | Array of read 2 FASTQ files representing a single GEX 10x library.| Array[File] |
| gex_i1_fastq | Optional array of index FASTQ files representing a single GEX 10x library; multiplexed samples are not currently supported, but the file may be passed to the pipeline. | Array[File] |
| tar_star_reference | TAR file containing a species-specific reference genome and GTF for Optimus (GEX) pipeline. | File | 
| ref_genome_fasta | Genome FASTA file used for building the indices. | File |
| mt_genes | Optional file for the Optimus (GEX) pipeline containing mitochondrial gene names used for metric calculation; default assumes 'mt' prefix in GTF (case insensitive). | File |
| counting_mode | Optional string that determines whether the Optimus (GEX) pipeline should be run in single-cell mode (sc_rna) or single-nucleus mode (sn_rna); default is "sn_rna". | String |
| tenx_chemistry_version | Optional integer for the Optimus (GEX) pipeline specifying the 10x version chemistry the data was generated with; validated by examination of the first read 1 FASTQ file read structure; default is "3". | Integer |
| emptydrops_lower | Optional threshold for UMIs for the Optimus (GEX) pipeline that empty drops tool should consider for determining cell; data below threshold is not removed; default is "100". | Integer |
| force_no_check | Optional boolean for the Optimus (GEX) pipeline indicating if the pipeline should perform checks; default is "false". | Boolean |
| ignore_r1_read_length | Optional boolean for the Optimus (GEX) pipeline indicating if the pipeline should ignore barcode chemistry check; if "true", the workflow will not ensure the `10x_chemistry_version` input matches the chemistry in the read 1 FASTQ; default is "false". | Boolean |
| star_strand_mode | Optional string for the Optimus (GEX) pipeline for performing STARsolo alignment on forward stranded, reverse stranded, or unstranded data; default is "Forward". | String |
| count_exons | Optional boolean for the Optimus (GEX) pipeline indicating if the workflow should calculate exon counts **when in single-nucleus (sn_rna) mode**; if "true" in sc_rna mode, the workflow will return an error; default is "false". | Boolean |
| gex_whitelist | Optional file containing the list of valid barcodes for 10x multiome GEX data; default is "gs://broad-gotc-test-storage/Multiome/input/737K-arc-v1_gex.txt". | File |
| atac_r1_fastq | Array of read 1 paired-end FASTQ files representing a single 10x multiome ATAC library. | Array[File] |
| atac_r2_fastq | Array of barcodes FASTQ files representing a single 10x multiome ATAC library. | Array[File] |
| atac_r3_fastq | Array of read 2 paired-end FASTQ files representing a single 10x multiome ATAC library. | Array[File] |
| tar_bwa_reference | TAR file containing the reference index files for BWA-mem alignment for the ATAC pipeline. | File | 
| chrom_sizes | File containing the genome chromosome sizes; used to calculate ATAC fragment file metrics. | File |
| adapter_seq_read1 | Optional string describing the adapter sequence for ATAC read 1 paired-end reads to be used during adapter trimming with Cutadapt; default is "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG". | String |
| adapter_seq_read3 | Optional string describing the adapter sequence for ATAC read 2 paired-end reads to be used during adapter trimming with Cutadapt; default is "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG". | String |
| atac_whitelist | Optional file containing the list of valid barcodes for 10x multiome ATAC adata; default is "gs://broad-gotc-test-storage/Multiome/input/737K-arc-v1_atac.txt". | File |


## Tasks

The Multiome workflow calls two subworkflows, which are described briefly in the table below. For more details on each subworkflow, including the tasks that they call, see the documentation linked in the table.

| Subworkflow | Software | Description | 
| ----------- | -------- | ----------- |
| ATAC ([WDL](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/multiome/atac.wdl) and [documentation](../ATAC/README)) | fastqprocess, bwa-mem, SnapATAC2 | Workflow used to analyze 10x single-cell ATAC data. |
| Optimus ([WDL](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/optimus/Optimus.wdl) and [documentation](../Optimus_Pipeline/README)) | fastqprocess, STARsolo, Emptydrops | Workflow used to analyze 10x single-cell GEX data. |
| JoinMultiomeBarcodes as JoinBarcodes ([WDL](https://github.com/broadinstitute/warp/blob/develop/tasks/skylab/H5adUtils.wdl)) | Python3 | Task that adds an extra column to the Optimus metrics `h5ad.obs` property that lists the respective ATAC barcodes for each gene expression barcode. It also adds an extra column to the ATAC metrics `h5ad.obs` property to link ATAC barcodes to gene expression barcodes. | 

## Outputs

| Output variable name | Filename, if applicable | Output format and description |
|--- | --- | --- | 
| multiome_pipeline_version_out | N.A. | String describing the version of the Multiome pipeline used. |
| bam_aligned_output_atac | `<input_id>_atac.bam` | BAM file containing aligned reads from ATAC workflow. |
| fragment_file_atac | `<input_id>_atac.fragments.tsv` | TSV file containing fragment start and stop coordinates per barcode. | 
| snap_metrics_atac | `<input_id>_atac.metrics.h5ad` | h5ad (Anndata) file containing per-barcode metrics from SnapATAC2. Also contains the equivalent gene expression barcode for each ATAC barcode in the `gex_barcodes` column of the `h5ad.obs` property. |
| genomic_reference_version_gex | `<reference_version>.txt` | File containing the Genome build, source and GTF annotation version. |
| bam_gex | `<input_id>_gex.bam` | BAM file containing aligned reads from Optimus workflow. |
| matrix_gex | `<input_id>_gex_sparse_counts.npz` | NPZ file containing raw gene by cell counts. |
| matrix_row_index_gex | `<input_id>_gex_sparse_counts_row_index.npy` | NPY file containing the row indices. |
| matrix_col_index_gex | `<input_id>_gex_sparse_counts_col_index.npy` | NPY file containing the column indices. |
| cell_metrics_gex | `<input_id>_gex.cell_metrics.csv.gz` | CSV file containing the per-cell (barcode) metrics. |
| gene_metrics_gex | `<input_id>_gex.gene_metrics.csv.gz` | CSV file containing the per-gene metrics. |
| cell_calls_gex | `<input_id>_gex.emptyDrops` | TSV file containing the EmptyDrops results when the Optimus workflow is run in sc_rna mode. |
| h5ad_output_file_gex | `<input_id>_gex.h5ad` | h5ad (Anndata) file containing the raw cell-by-gene count matrix, gene metrics, cell metrics, and global attributes. Also contains equivalent ATAC barcode for each gene expression barcode in the `atac_barcodes` column of the `h5ad.obs` property. |

## Versioning and testing

All Multiome pipeline releases are documented in the [Multiome changelog](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/multiome/Multiome.changelog.md) and tested using [plumbing and scientific test data](https://github.com/broadinstitute/warp/tree/master/pipelines/skylab/multiome/test_inputs). To learn more about WARP pipeline testing, see [Testing Pipelines](https://broadinstitute.github.io/warp/docs/About_WARP/TestingPipelines).

## Citing the Multiome Pipeline
Please identify the pipeline in your methods section using the Multiome Pipeline's [SciCrunch resource identifier](https://scicrunch.org/resources/data/record/nlx_144509-1/SCR_024217/resolver?q=SCR_024217&l=SCR_024217&i=rrid:scr_024217).
* Ex: *Multiome Pipeline (RRID:SCR_024217)*

## Consortia support
This pipeline is supported by the [BRAIN Initiative](https://braininitiative.nih.gov/) (BICCN and BICAN). 

If your organization also uses this pipeline, we would like to list you! Please reach out to us by contacting the [WARP Pipeline Development team](mailto:warp-pipelines-help@broadinstitute.org).

## Acknowledgements
We are immensely grateful to the members of the BRAIN Initiative (BICAN Sequencing Working Group) and SCORCH for their invaluable and exceptional contributions to this pipeline. Our heartfelt appreciation goes to Alex Dobin, Aparna Bhaduri, Alec Wysoker, Anish Chakka, Brian Herb, Daofeng Li, Fenna Krienen, Guo-Long Zuo, Jeff Goldy, Kai Zhang, Khalid Shakir, Bo Li, Mariano Gabitto, Michael DeBerardine, Mengyi Song, Melissa Goldman, Nelson Johansen, James Nemesh, and Theresa Hodges for their unwavering dedication and remarkable efforts. 

## Feedback

Please help us make our tools better by contacting the [WARP Pipeline Development team](mailto:warp-pipelines-help@broadinstitute.org) for pipeline-related suggestions or questions.
