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

The workflow is a wrapper WDL script that consists of two subworkflows: the Optimus workflow for GEX data and the ATAC workflow for single-cell ATAC data. The GEX component performs barcode and UMI correction, aligns reads to the genome, and produces both quality metrics per barcode and gene and a raw cell-by-gene count matrix. The ATAC component corrects cell barcodes, aligns reads to the genome, and producesa fragment file as well as per barcode metrics. 

Read more about the GEX workflow in the Optimus overview for the ATAC workflow in the ATAC workflow. 

The wrapper WDL is available in the WARP repository (see the [code here](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/multiome/Multiome.wdl)).

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

To download the latest Multiome release, see the release tags prefixed with "Multiome" on the WARP [releases page](https://github.com/broadinstitute/warp/releases). All Multiome pipeline releases are documented in the [Multiome changelog](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/multiome/Multiome.changelog.md). 

To discover and search releases, use the WARP command-line tool [Wreleaser](https://github.com/broadinstitute/warp/tree/master/wreleaser).

If you’re running an ultiomeus workflow version prior to the latest release, the accompanying documentation for that release may be downloaded with the source code on the WARP [releases page](https://github.com/broadinstitute/warp/releases) (see the source code folder).

Multiome can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH compliant, flexible workflow management system that supports multiple computing platforms. The workflow can also be run in [Terra](https://app.terra.bio), a cloud-based analysis platform. 
## Input variables

| Variable name | Description | Optional attributes (when applicable) |
| --- | --- | --- |
| counting_mode | Optional string that demarcates whether the pipeline should be run in single-cell mode (sc_rna) or single-nucleus mode (sn_rna) | "sn_rna" is set by default. |
| gex_r1_fastq | Array of read 1 FASTQ files representing a single 10x library. |
| gex_r2_fastq | Array of read 2 FASTQ files representing a single 10x library.| N.A. |
| gex_i1_fastq | Array of index FASTQ files representing a single 10x library. The pipeline does not currently support multiplexed samples, but the file can be optionally passed through. | N.A. |
| input_id | String demarcating the sampe ID or name; used to name the GEX output files. | N.A. |
| output_bam_basename | String demarcating the sample ID or name used to name the output BAM file | Set to input_id by default. |
| tar_star_reference | Tar file containing reference index files for STARsolo alignment. | N.A. |
| annotations_gtf | A GTF file with gene annotations to use for GEX cell metric calculation and ATAC fragment metrics. This GTF should match the GTF used to build the STAR aligner. | N.A. |
| ref_genome_fasta | Genome FASTA file used for building the indices. | N.A. |
| mt_genes | Optional file that lists mitochondria genes to be used for cell metric calculation. | Default assumes 'mt' prefix in GTF (case insensitive). |
| tenx_chemistry_version | Integer specifying the 10x barcode read structure. Multiome uses same barcode readstructure as 10x v3. | Set to 3 by default because 10x Multiome data uses same barcode read structure as v3. |
| emptydrops_lower | A threshold for number of UMIs empty drops tool should consider for determining cell. Returns a 'yes' or 'no' value, but does not remove data.  | Default is set to 100. |
| force_no_check | Boolean indicating if the pipeline should performs checks. | Default is set to false |
| ignore_r1_read_length | Boolean to indicate whether to ignore checks for barcode chemistry | Default set to false |
| star_strand_mode | String indicating which strand option to use for STARsolo alginment. This should match the parameters specified in STAR documentation. | Default is set "Forward" to match 10x chemistry|
| count_exons | Boolean to indicate whether to run STAR in both gene_full and gene modes. | Default is set to false. |
| gex_whitelist | File containing list of valid barcodes for 10x multiome gene expression data. | "gs://broad-gotc-test-storage/Multiome/input/737K-arc-v1_gex.txt" |
| atac_r1_fastq | Array[File] | N.A. |
| atac_r2_fastq | Array[File] | N.A. |
| atac_r3_fastq | Array[File] | N.A. |
| output_base_name | String | N.A. |
| tar_bwa_reference | File | N.A. |
| chrom_sizes | File | N.A. |
| adapter_seq_read1 | String | "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG" |
| adapter_seq_read3 | String | "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG" |
| atac_whitelist | File containing list of valid barcodes for 10x multiome ATAC adata. | "gs://broad-gotc-test-storage/Multiome/input/737K-arc-v1_atac.txt" |


## Tasks
| Task name and WDL link | Tool | Software | Description | 
| --- | --- | --- | --- | 
| 

## Output Variables

| Variable Name         | Description                                    |
|----------------------|------------------------------------------------|
| bam_aligned_output    | Aligned BAM file                               |
| fragment_file         | Fragment file                                  |
| snap_metrics          | Snap metrics file                              |