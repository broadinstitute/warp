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
| Assay type | single-nucleus methylome and chromatin contact (snM3C) sequencing data | [Lee et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6765423/) |
| Overall workflow | Read alignment and chromatin contact calling |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Genomic Reference Sequence | GRCh38 human genome primary sequence | GENCODE [human reference files](https://www.gencodegenes.org/human/release_43.html)|
| Aligner | HISAT-3N | [Zhang at al. 2021](https://genome.cshlp.org/content/31/7/1290) |
| Data input file format | File format in which sequencing data is provided | [FASTQ](https://academic.oup.com/nar/article/38/6/1767/3112533) |
| Data output file format | File formats in which snM3C output is provided | TSV, [FASTQ](https://academic.oup.com/nar/article/38/6/1767/3112533), [BAM](http://samtools.github.io/hts-specs/), and [ALLC](https://lhqing.github.io/ALLCools/intro.html) |


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

1. Demultiplexes, sorts, and trims reads.
2. Aligns paired-end reads.
3. Separates unmapped, uniquely aligned, multi-aligned reads.
4. Splits unmapped reads by enzyme cut sites.
5. Aligns unmapped, single-end reads.
6. Removes overlapping reads.
7. Merges mapped reads from single- and paired-end alignments.
8. Calls chromatin contacts.
9. Removes duplicate reads.
10. Creates ALLC file.
11. Creates summary output file.

The tools each snM3C task employs are detailed in the table below. 

To see specific tool parameters, select the [workflow WDL link](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/snM3C/snM3C.wdl); then find the task and view the `command {}` section of the task in the WDL script. To view or use the exact tool software, see the task's Docker image which is specified in the task WDL `# runtime values` section as `docker: `. More details about these tools and parameters can be found in the [YAP documentation](https://hq-1.gitbook.io/mc/).

| Task name | Tool | Software | Description |
| --- | --- | --- | --- |
| Demultiplexing | Cutadapt | [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) | Performs demultiplexing to cell-level FASTQ files based on random primer indices. |
| Sort_and_trim_r1_and_r2 | Cutadapt | [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) | Sorts, filters, and trims reads using the `r1_adapter`, `r2_adapter`, `r1_left_cut`, `r1_right_cut`, `r2_left_cut`, and `r2_right_cut` input parameters. |
| Hisat_3n_pair_end_mapping_dna_mode | HISAT-3N | [HISAT-3N](https://daehwankimlab.github.io/hisat2/hisat-3n/) | Performs paired-end read alignment. |
| Separate_unmapped_reads | [hisat3n_general.py](https://github.com/lhqing/cemba_data/blob/788e83cd66f3b556bdfacf3485bed9500d381f23/cemba_data/hisat3n/hisat3n_general.py) | python3 | Imports a custom python3 script developed by Hanqing Liu and calls the `separate_unique_and_multi_align_reads()` function to separate unmapped, uniquely aligned, multi-aligned reads from HISAT-3N BAM file; unmapped reads are stored in an unmapped FASTQ file and uniquely and multi-aligned reads are stored in separate BAM files. |
| Split_unmapped_reads | [hisat3n_m3c.py](https://github.com/lhqing/cemba_data/blob/bf6248239074d0423d45a67d83da99250a43e50c/cemba_data/hisat3n/hisat3n_m3c.py) | python3 | Imports a custom python3 script developed by Hanqing Liu and calls the `split_hisat3n_unmapped_reads()` function to split the unmapped reads FASTQ file by all possible enzyme cut sites and output new R1 and R2 FASTQ files. |
| Hisat_single_end_r1_r2_mapping_dna_mode_and_merge_sort_split_reads_by_name | HISAT-3N | [HISAT-3N](https://daehwankimlab.github.io/hisat2/hisat-3n/) | Performs single-end alignment of unmapped reads to maximize read mapping. |
| remove_overlap_read_parts | [hisat3n_m3c.py](https://github.com/lhqing/cemba_data/blob/bf6248239074d0423d45a67d83da99250a43e50c/cemba_data/hisat3n/hisat3n_m3c.py) | python3 | Imports a custom python3 script developed by Hanqing Liu and calls the `remove_overlap_read_parts()` function to remove overlapping reads from the split alignment BAM file produced during single-end alignment. |
| merge_original_and_split_bam_and_sort_all_reads_by_name_and_position | merge, sort | [samtools](https://www.htslib.org/) | Merges and sorts all mapped reads from the paired-end and single-end alignments; creates a position-sorted BAM file and a name-sorted BAM file. |
| call_chromatin_contacts | [hisat3n_m3c.py](https://github.com/lhqing/cemba_data/blob/bf6248239074d0423d45a67d83da99250a43e50c/cemba_data/hisat3n/hisat3n_m3c.py) | python3 | Imports a custom python3 script developed by Hanqing Liu and calls the `call_chromatin_contacts()` function to call chromatin contacts from the name-sorted, merged BAM file; reads are considered chromatin contacts if they are greater than 2,500 base pairs apart. |
| dedup_unique_bam_and_index_unique_bam | MarkDuplicates | [Picard](https://broadinstitute.github.io/picard/) | Removes duplicate reads from the position-sorted, merged BAM file. |
| unique_reads_allc | bam-to-allc | [ALLCools](https://lhqing.github.io/ALLCools/intro.html) | Creates an ALLC file with a list of methylation points. |
| unique_reads_cgn_extraction | extract-allc | [ALLCools](https://lhqing.github.io/ALLCools/intro.html) | Creates an ALLC file containing methylation contexts. |
| summary | [summary.py](https://github.com/lhqing/cemba_data/blob/788e83cd66f3b556bdfacf3485bed9500d381f23/cemba_data/hisat3n/summary.py) | python3 | Imports a custom python3 script developed by Hanqing Liu and calls the `snm3c_summary()` function to generate a single, summary file for the pipeline in TSV format; contains trimming, mapping, deduplication, chromatin contact, and AllC site statistics. |

#### 1. Demultiplexes, sorts, and trims reads
In the first step of the pipeline (`Demultiplexing`), raw sequencing reads are demultiplexed by random primer index into cell-level FASTQ files using [Cutadapt](https://cutadapt.readthedocs.io/en/stable/). For more information on barcoding, see the [YAP documentation](https://hq-1.gitbook.io/mc/tech-background/barcoding#two-round-of-barcoding). 

After demultiplexing, the pipeline uses [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) to sort, filter, and trim reads in the `Sort_and_trim_r1_and_r2` task. The R1 and R2 adapter sequences are removed, along with the number of bases specified by the `r1_left_cut`, `r1_right_cut`, `r2_left_cut`, and `r2_right_cut` input parameters. Any reads shorter than the specified `min_read_length` are filtered out in this step.

#### 2. Aligns paired-end reads
In the next step of the pipeline, the `Hisat_3n_pair_end_mapping_dna_mode` task uses [HISAT-3N](https://daehwankimlab.github.io/hisat2/hisat-3n/) to perform paired-end read alignment to a reference genome FASTA file (`genome_fa`) and outputs an aligned BAM file. Additionally, the task outputs a stats file and a text file containing the genomic reference version used.

#### 3. Separates unmapped, uniquely aligned, multi-aligned reads
After paired-end alignment, the pipeline calls the `Separate_unmapped_reads` task, which imports a custom python3 script ([hisat3n_general.py](https://github.com/lhqing/cemba_data/blob/788e83cd66f3b556bdfacf3485bed9500d381f23/cemba_data/hisat3n/hisat3n_general.py)) developed by Hanqing Liu. The task calls the script's `separate_unique_and_multi_align_reads()` function to separate unmapped, uniquely aligned, and multi-aligned reads from the HISAT-3N BAM file. Three new files are output from this step of the pipeline: 

1. A FASTQ file that contains the unmapped reads (`unmapped_fastq_tar`)
2. A BAM file that contains the uniquely aligned reads (`unique_bam_tar`)
3. A BAM file that contains the multi-aligned reads (`multi_bam_tar`)

#### 4. Splits unmapped reads by enzyme cut sites
The `Split_unmapped_reads` task imports a custom python3 script ([hisat3n_m3c.py](https://github.com/lhqing/cemba_data/blob/bf6248239074d0423d45a67d83da99250a43e50c/cemba_data/hisat3n/hisat3n_m3c.py)) developed by Hanqing Liu and calls the script's `split_hisat3n_unmapped_reads()` function. This splits the FASTQ file containing the unmapped reads by all possible enzyme cut sites and outputs new R1 and R2 files. 

#### 5. Aligns unmapped, single-end reads
In the next step of the pipeline, the `Hisat_single_end_r1_r2_mapping_dna_mode_and_merge_sort_split_reads_by_name ` task uses [HISAT-3N](https://daehwankimlab.github.io/hisat2/hisat-3n/) to perform single-end read alignment of the previously unmapped reads to maximize read mapping and outputs a single, aligned BAM file.

#### 6. Removes overlapping reads
After the second alignment step, the pipeline calls the `remove_overlap_read_parts ` task, which imports a custom python3 script ([hisat3n_m3c.py](https://github.com/lhqing/cemba_data/blob/bf6248239074d0423d45a67d83da99250a43e50c/cemba_data/hisat3n/hisat3n_m3c.py)) developed by Hanqing Liu. The task calls the script's `remove_overlap_read_parts()` function to remove overlapping reads from the BAM file produced during single-end alignment and output another BAM file.

#### 7. Merges mapped reads from single- and paired-end alignments
The `merge_original_and_split_bam_and_sort_all_reads_by_name_and_position` task uses [samtools](https://www.htslib.org/) to merge and sort all of the mapped reads from the paired-end and single-end alignments into a single BAM file. The BAM file is output as both a position-sorted and a name-sorted BAM file.

#### 8. Calls chromatin contacts
In the `call_chromatin_contacts` task, the pipeline imports a custom python3 script ([hisat3n_m3c.py](https://github.com/lhqing/cemba_data/blob/bf6248239074d0423d45a67d83da99250a43e50c/cemba_data/hisat3n/hisat3n_m3c.py)) developed by Hanqing Liu. The task calls the script's `call_chromatin_contacts()` function to call chromatin contacts from the name-sorted, merged BAM file. If reads are greater than 2,500 base pairs apart, they are considered chromatin contacts. If reads are less than 2,500 base pairs apart, they are considered the same fragment. 

#### 9. Removes duplicate reads
After calling chromatin contacts, the `dedup_unique_bam_and_index_unique_bam` task uses Picard's MarkDuplicates tool to remove duplicate reads from the position-sorted, merged BAM file and output a deduplicated BAM file.

#### 10. Creates ALLC file
The `unique_reads_allc` task uses the [ALLCools](https://lhqing.github.io/ALLCools/intro.html) `bam-to-allc` function to create an ALLC file from the deduplicated BAM file that contains a list of methylation points. The `num_upstr_bases` and `num_downstr_bases` input parameters are used to define the number of bases upstream and downstream of the C base to include in the ALLC context column.

Next, the `unique_reads_cgn_extraction` task uses the [ALLCools](https://lhqing.github.io/ALLCools/intro.html) `extract-allc` function to extract methylation contexts from the input ALLC file and output a second ALLC file that can be used to generate an [MCDS file](https://github.com/lhqing/allcools_doc/blob/master/tech-background/file-formats.md#mcds-file). 

#### 11. Creates summary output file
In the last step of the pipeline, the `summary` task imports a custom python3 script ([summary.py](https://github.com/lhqing/cemba_data/blob/788e83cd66f3b556bdfacf3485bed9500d381f23/cemba_data/hisat3n/summary.py)) developed by Hanqing Liu. The task calls the script's `snm3c_summary()` function to generate a single, summary file for the pipeline in TSV format; contains trimming, mapping, deduplication, chromatin contact, and AllC site statistics. This is the main output of the pipeline.

## Outputs

The following table lists the output variables and files produced by the pipeline.

| Output name | Filename, if applicable | Output format and description |
| ------ | ------ | ------ |
| MappingSummary | `<plate_id>_MappingSummary.csv.gz` | Mapping summary file in CSV format. |

<!--- describe outputs --->

        

| allcFiles | Tarred file containing allc files |
| allc_CGNFiles| Tarred file containing CGN context-specific allc files | 
| bamFiles | Tarred file containing cell-level aligned BAM files |
| detail_statsFiles | Tarred file containing detail stats files | 
| hicFiles | Tarred file containing Hi-C files |

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


## Versioning

All snM3C pipeline releases are documented in the [pipeline changelog](https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/snM3C/snM3C.changelog.md).

## Consortia support
This pipeline is supported by the [BRAIN Initiative](https://braininitiative.nih.gov/) (BICCN and BICAN). 

If your organization also uses this pipeline, we would like to list you! Please reach out to us by contacting the [WARP Pipeline Development team](mailto:warp-pipelines-help@broadinstitute.org).

## Feedback

For questions, suggestions, or feedback related to the snM3C pipeline, please contact [the WARP team](mailto:warp-pipelines-help@broadinstitute.org). Your feedback is valuable for improving the pipeline and addressing any issues that may arise during its usage.