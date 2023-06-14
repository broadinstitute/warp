---
sidebar_position: 1
---

# Multiome Overview

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [Multiome v1.0.1](https://github.com/broadinstitute/warp/releases) | June, 2023 | Kaylee Mathews | Please file GitHub issues in warp or contact the [WARP Pipeline Development team](mailto:warp-pipelines-help@broadinstitute.org) |

![Multiome_diagram]()

## Introduction to the Multiome workflow

Multiome is an open-source, cloud-optimized pipeline developed in collaboration with members of the [BRAIN Initiative](https://braininitiative.nih.gov/) (BICCN and BICAN), including the [Allen Institute for Brain Science](https://alleninstitute.org/division/brain-science/), [Neuroscience MultiOmic Archive](https://nemoarchive.org/), Kai Zhang ([SnapATAC2](https://kzhang.org/SnapATAC2/index.html)), and Alex Dobin ([STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md)). It supports the processing of 10x 3' single-cell and single-nucleus gene expression (GEX) and chromatin accessibility (ATAC) data generated with the [10x Genomics Multiome assay](https://www.10xgenomics.com/products/single-cell-multiome-atac-plus-gene-expression).

The workflow is a wrapper WDL script that calls two subworkflows: the [Optimus workflow](../Optimus_Pipeline/README) for GEX data and the [ATAC workflow](../ATAC/README) for single-cell ATAC data. 

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
| counting_mode | Optional string that determines whether the pipeline should be run in single-cell mode (sc_rna) or single-nucleus mode (sn_rna); default is "sn_rna". | String |
| gex_r1_fastq | Array of read 1 FASTQ files representing a single 10x library. | Array[File] |
| gex_r2_fastq | Array of read 2 FASTQ files representing a single 10x library.| Array[File] |
| gex_i1_fastq | Optional array of index FASTQ files representing a single 10x library; multiplexed samples are not currently supported, but the file may be passed to the pipeline. | Array[File] |
| input_id | Unique identifier describing the biological sample or replicate that corresponds with the FASTQ files; can be a human-readable name or UUID; used to name the GEX output files. | String |
| output_bam_basename | Used as basename for output BAM file; default is `input_id`. | String |
| tar_star_reference | TAR file containing a species-specific reference genome and GTF. | File | 
| annotations_gtf | GTF file containing gene annotations used for GEX cell metric calculation and ATAC fragment metrics; must match the GTF used to build the STAR aligner. | File |
| ref_genome_fasta | Genome FASTA file used for building the indices. | File |
| mt_genes | Optional file containing mitochondrial gene names used for metric calculation; default assumes 'mt' prefix in GTF (case insensitive). | File |
| tenx_chemistry_version | Optional integer specifying the 10x version chemistry the data was generated with; validated by examination of the first read 1 FASTQ file read structure; default is "3". | Integer |
| emptydrops_lower | Optional threshold for UMIs that empty drops tool should consider for determining cell; data below threshold is not removed; default is "100". | Integer |
| force_no_check | Optional boolean indicating if the pipeline should perform checks; default is "false". | Boolean |
| ignore_r1_read_length | Optional boolean indicating if the pipeline should ignore barcode chemistry check; if "true", the workflow will not ensure the `10x_chemistry_version` input matches the chemistry in the read 1 FASTQ; default is "false". | Boolean |
| star_strand_mode | Optional string for performing STARsolo alignment on forward stranded, reverse stranded, or unstranded data; default is "Forward". | String |
| count_exons | Optional boolean indicating if the workflow should calculate exon counts **when in single-nucleus (sn_rna) mode**; if "true" in sc_rna mode, the workflow will return an error; default is "false". | Boolean |
| gex_whitelist | Optional file containing the list of valid barcodes for 10x multiome GEX data; default is "gs://broad-gotc-test-storage/Multiome/input/737K-arc-v1_gex.txt". | File |
| atac_r1_fastq | Array of read 1 paired-end FASTQ files representing a single 10x multiome ATAC library. | Array[File] |
| atac_r2_fastq | Array of barcodes FASTQ files representing a single 10x multiome ATAC library. | Array[File] |
| atac_r3_fastq | Array of read 2 paired-end FASTQ files representing a single 10x multiome ATAC library. | Array[File] |
| output_base_name | Used as basename for output ATAC files. | String |
| tar_bwa_reference | TAR file containing the reference index files for BWA-mem alignment. | File | 
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

## Outputs

| Output variable name | Filename, if applicable | Output format and description |
|--- | --- | --- | 
| bam_aligned_output | `<input_id>.bam` | BAM file containing aligned reads from ATAC workflow. |
| fragment_file | `<input_id>.fragments.tsv` | TSV file containing fragment start and stop coordinates per barcode. | 
| snap_metrics | `<input_id>.metrics.h5ad` | h5ad (Anndata) file containing per-barcode metrics from SnapATAC2. |
| pipeline_version_out | N.A. | String describing the Optimus pipeline version used. |
| genomic_reference_version | `<reference_version>.txt` | File containing the Genome build, source and GTF annotation version. |
| bam | `<input_id>.bam` | BAM file containing aligned reads from Optimus workflow. |
| matrix | `<input_id>_sparse_counts.npz` | NPZ file containing raw gene by cell counts. |
| matrix_row_index | `<input_id>_sparse_counts_row_index.npy` | NPY file containing the row indices. |
| matrix_col_index | `<input_id>_sparse_counts_col_index.npy` | NPY file containing the column indices. |
| cell_metrics | `<input_id>.cell_metrics.csv.gz` | CSV file containing the per-cell (barcode) metrics. |
| gene_metrics | `<input_id>.gene_metrics.csv.gz` | CSV file containing the per-gene metrics. |
| cell_calls | `<input_id>.emptyDrops` | TSV file containing the EmptyDrops results when the Optimus workflow is run in sc_rna mode. |
| h5ad_output_file | `<input_id>.h5ad` | h5ad (Anndata) file containing the raw cell-by-gene count matrix, gene metrics, cell metrics, and global attributes. |

## Versioning and testing

All Multiome pipeline releases are documented in the [Multiome changelog](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/multiome/Multiome.changelog.md) and tested using [plumbing and scientific test data](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/multiome/test_inputs/test_data_overview.md). To learn more about WARP pipeline testing, see [Testing Pipelines](https://broadinstitute.github.io/warp/docs/About_WARP/TestingPipelines).

## Consortia support
This pipeline is supported by the [BRAIN Initiative](https://braininitiative.nih.gov/) (BICCN and BICAN). 

If your organization also uses this pipeline, we would like to list you! Please reach out to us by contacting the [WARP Pipeline Development team](mailto:warp-pipelines-help@broadinstitute.org).

## Feedback

Please help us make our tools better by contacting the [WARP Pipeline Development team](mailto:warp-pipelines-help@broadinstitute.org) for pipeline-related suggestions or questions.