---
sidebar_position: 1
slug: /Pipelines/snM3C/README
---
# Single Nucleus Methyl-Seq and Chromatin Capture (snM3C) Overview

| Pipeline Version | Date Updated | Documentation Authors | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [snM3C_v1.0.0](https://github.com/broadinstitute/warp/releases) | January, 2024 | [Kaylee Mathews](mailto:warp-pipelines-help@broadinsitute.org) | Please file GitHub issues in the [WARP repository](https://github.com/broadinstitute/warp/issues) |


## Introduction to snM3C

The Single Nucleus Methly-Seq and Chromatin Capture (snM3C) workflow is an open-source, cloud-optimized computational workflow for processing single-nucleus methylome and chromatin contact (snM3C) sequencing data. The workflow is designed to demultiplex and align raw sequencing reads, call chromatin contacts, and generate summary metrics. 

The workflow is developed in collaboration Hanqing Liu and the laboratory of Joseph Ecker. For more information about the snM3C tools and analysis, please see the [YAP documentation](https://hq-1.gitbook.io/mc/) or the [cemba_data](https://github.com/lhqing/cemba_data) GitHub repository created by Hanqing Liu.

## Quickstart table
The following table provides a quick glance at the Multiome pipeline features:

| Pipeline features | Description | Source |
|--- | --- | --- |
| Assay type | single-nucleus methylome and chromatin contact (snM3C) sequencing data | [Lee et al. 2020](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6765423/) |
<!---
| Overall workflow  | Barcode correction, read alignment, gene and fragment quanitification |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Genomic Reference Sequence | GRCh38 human genome primary sequence | GENCODE [human reference files](https://www.gencodegenes.org/human/release_43.html)|
| Gene annotation reference (GTF) | Reference containing gene annotations | GENCODE [human GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz) |
| Aligners | STARsolo (GEX), BWA-mem2 (ATAC) | [Kaminow et al. 2021](https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1), [Vasimuddin et al. 2019](https://ieeexplore.ieee.org/document/8820962) |
| Transcript and fragment quantification | STARsolo (GEX), SnapATAC2 (ATAC) | [Kaminow et al. 2021](https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1), [SnapATAC2](https://kzhang.org/SnapATAC2/) |
| Data input file format | File format in which sequencing data is provided | [FASTQ](https://academic.oup.com/nar/article/38/6/1767/3112533) |
| Data output file format | File formats in which Multiome output is provided | [BAM](http://samtools.github.io/hts-specs/) and [h5ad](https://anndata.readthedocs.io/en/latest/) |
--->


## Set-up

### snM3C installation

To download the latest snM3C release, see the release tags prefixed with "snM3C" on the WARP [releases page](https://github.com/broadinstitute/warp/releases). All snM3C pipeline releases are documented in the [snM3C changelog](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/snM3C/snM3C.changelog.md). 

To discover and search releases, use the WARP command-line tool [Wreleaser](https://github.com/broadinstitute/warp/tree/master/wreleaser).

If you’re running a version of the snM3C workflow prior to the latest release, the accompanying documentation for that release may be downloaded with the source code on the WARP [releases page](https://github.com/broadinstitute/warp/releases) (see the source code folder `website/docs/Pipelines/snM3C`).

The snM3C workflow can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant, flexible workflow management system that supports multiple computing platforms. The workflow can also be run in [Terra](https://app.terra.bio), a cloud-based analysis platform. 


### Inputs

The snM3C workflow requires a JSON configuration file specifying the input files and parameters for the analysis. Example configuration files can be found in the snM3C [`test_inputs`](https://github.com/broadinstitute/warp/tree/develop/pipelines/skylab/snM3C/test_inputs) directory in the WARP repository.

#### Input descriptions

| Parameter | Description |
| ---| --- |
| fastq_input_read1 | Array of multiplexed FASTQ files for read 1. |
| fastq_input_read2 | Array of multiplexed FASTQ files for read 2. |
| random_primer_indexes | File containing random primer indexes. |
| plate_id | String specifying the plate ID. |
| tarred_index_files | File containing tarred index files for hisat-3 mapping. |
| genome_fa | File containing the reference genome in FASTA format. | 
| chromosome_sizes | File containing the genome chromosome sizes. |
| r1_adapter | Optional string describing the adapter sequence for read 1 paired-end reads to be used during adapter trimming with Cutadapt; default is "AGATCGGAAGAGCACACGTCTGAAC". |
| r2_adapter | Optional string describing the adapter sequence for read 2 paired-end reads to be used during adapter trimming with Cutadapt; default is  "AGATCGGAAGAGCGTCGTGTAGGGA". |
| r1_left_cut | Optional integer describing the number of bases to be trimmed from the beginning of read 1 with Cutadapt; default is 10. |
| r1_right_cut | Optional integer describing the number of bases to be trimmed from the end of read 1 with Cutadapt; default is 10. |
| r2_left_cut | Optional integer describing the number of bases to be trimmed from the beginning of read 2 with Cutadapt; default is 10. |
| r2_right_cut | Optional integer describing the number of bases to be trimmed from the end of read 2 with Cutadapt; default is 10. |
| min_read_length | Optional integer; if a read length is smaller than `min_read_length`, both paired-end reads will be discarded; default is 30.  |
| num_upstr_bases | Optional integer describing the number of bases upstream of the C base to include in ALLC file context column created using ALLCools; default is 0. |
| num_downstr_bases | Optional integer describing the number of bases upstream of the C base to include in ALLC file context column created using ALLCools; default is 2. |
| compress_level | Optional integer describing the compression level for the output ALLC file; default is 5. |


## snM3C tasks and tools
The workflow contains several tasks described below.

Overall, the snM3C workflow:
<!--- describe what it does
1.
--->

The tools each snM3C task employs are detailed in the table below. 

To see specific tool parameters, select the [workflow WDL link](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/snM3C/snM3C.wdl); then find the task and view the `command {}` section of the task in the WDL script. To view or use the exact tool software, see the task's Docker image which is specified in the task WDL `# runtime values` section as `docker: `. More details about these tools and parameters can be found in the [YAP documentation](https://hq-1.gitbook.io/mc/).

| Task name | Tool | Software | Description |
| --- | --- | --- | --- |
| Demultiplexing | Cutadapt | [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) | Performs demultiplexing to cell-level FASTQ files. |
| Sort_and_trim_r1_and_r2 | Cutadapt | [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) | Sorts, filters, and trims reads. |
| Hisat_3n_pair_end_mapping_dna_mode | HISAT-3N | [HISAT-3N](https://daehwankimlab.github.io/hisat2/hisat-3n/) | Performs paired-end read alignment. |
<!--- need more detail
| Separate_unmapped_reads | tool | software | description |
| Split_unmapped_reads | tool | software | description |
| Hisat_single_end_r1_r2_mapping_dna_mode_and_merge_sort_split_reads_by_name | HISAT-3N | [HISAT-3N](https://daehwankimlab.github.io/hisat2/hisat-3n/) | Performs paired-end read alignment. |
| remove_overlap_read_parts | python3 | python3 | description |
| merge_original_and_split_bam_and_sort_all_reads_by_name_and_position | tool | software | description |
| call_chromatin_contacts | tool | software | description |
| dedup_unique_bam_and_index_unique_bam | tool | software | description |
| unique_reads_allc | bam-to-allc | ALLCools | description |
| unique_reads_cgn_extraction | extract-allc | ALLCools | description |
| summary | tool | software | description |

| Mapping | hisat-3 | hisat-3 | Performs trimming, alignment and calling chromatin contacts with a [custom snakemake](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/snM3C/Config%20files/Snakemake-file/Snakefile) file developed by Hanqing Liu. |
--->

<!--- describe tasks --->

## Outputs

The following table lists the output variables and files produced by the pipeline.

| Output name | Filename, if applicable | Output format and description |
| ------ | ------ | ------ |
<!--- add detail/update 
        File MappingSummary = summary.mapping_summary
        File trimmed_stats = Sort_and_trim_r1_and_r2.trim_stats_tar
        File r1_trimmed_fq = Sort_and_trim_r1_and_r2.r1_trimmed_fq_tar
        File r2_trimmed_fq = Sort_and_trim_r1_and_r2.r2_trimmed_fq_tar
        File hisat3n_stats_tar = Hisat_3n_pair_end_mapping_dna_mode.hisat3n_paired_end_stats_tar
        File hisat3n_bam_tar = Hisat_3n_pair_end_mapping_dna_mode.hisat3n_paired_end_bam_tar
        File unique_bam_tar = Separate_unmapped_reads.unique_bam_tar
        File multi_bam_tar = Separate_unmapped_reads.multi_bam_tar
        File unmapped_fastq_tar = Separate_unmapped_reads.unmapped_fastq_tar
        File split_fq_tar = Split_unmapped_reads.split_fq_tar
        File merge_sorted_bam_tar = Hisat_single_end_r1_r2_mapping_dna_mode_and_merge_sort_split_reads_by_name.merge_sorted_bam_tar
        File name_sorted_bams = merge_original_and_split_bam_and_sort_all_reads_by_name_and_position.name_sorted_bam
        File pos_sorted_bams = merge_original_and_split_bam_and_sort_all_reads_by_name_and_position.position_sorted_bam
        File remove_overlap_read_parts_bam_tar = remove_overlap_read_parts.output_bam_tar
        File dedup_unique_bam_and_index_unique_bam_tar = dedup_unique_bam_and_index_unique_bam.output_tar
        File unique_reads_cgn_extraction_allc = unique_reads_cgn_extraction.output_allc_tar
        File unique_reads_cgn_extraction_tbi = unique_reads_cgn_extraction.output_tbi_tar
        File chromatin_contact_stats = call_chromatin_contacts.chromatin_contact_stats
        File reference_version = Hisat_3n_pair_end_mapping_dna_mode.reference_version
        
| mappingSummary | Mapping summary file in CSV format |
| allcFiles | Tarred file containing allc files |
| allc_CGNFiles| Tarred file containing CGN context-specific allc files | 
| bamFiles | Tarred file containing cell-level aligned BAM files |
| detail_statsFiles | Tarred file containing detail stats files | 
| hicFiles | Tarred file containing Hi-C files |
--->


## Versioning

All snM3C pipeline releases are documented in the [pipeline changelog](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/snM3C/snM3C.changelog.md).

## Consortia support
This pipeline is supported by the [BRAIN Initiative](https://braininitiative.nih.gov/) (BICCN and BICAN). 

If your organization also uses this pipeline, we would like to list you! Please reach out to us by contacting the [WARP Pipeline Development team](mailto:warp-pipelines-help@broadinstitute.org).

## Feedback

For questions, suggestions, or feedback related to the snM3C pipeline, please contact [the WARP team](mailto:warp-pipelines-help@broadinstitute.org). Your feedback is valuable for improving the pipeline and addressing any issues that may arise during its usage.