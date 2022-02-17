---
sidebar_position: 2
---

# RNA with UMIs v1.0.1 Methods

Below we provide an example methods section for publications using the RNA with UMIs pipeline. For the complete pipeline documentation, see the [RNA with UMIs Overview](./README.md).

## Methods

Data preprocessing, gene counting, and metric calculation were performed using the RNA with UMIs v1.0.1 pipeline, which uses Picard v2.26.6, fgbio v1.4.0, STAR v2.7.10a, Samtools 1.11, UMI-tools v1.1.1, RNA-SeQC v2.4.2 with default tool parameters unless otherwise specified. Reference files are publicly available in the [Broad References](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references;tab=objects?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false) Google Bucket and are also listed in [example configuration files](https://github.com/broadinstitute/warp/tree/develop/pipelines/broad/rna_seq/test_inputs) in the in the WARP repository.

Paired-end FASTQ files were first converted to an unmapped BAM using Picard's FastqToSam tool with SORT_ORDER = unsorted. (If a read group unmapped BAM file is used as input for the pipeline, this step is skipped.) Unique molecular identifiers (UMIs) were extracted from the unmapped BAM using fgbio's ExtractUmisFromBam and stored in the RX read tag.

After the extraction of UMIs, the unmapped BAM file was aligned to the GENCODE GRCh38 (hg38) v34 (or GENCODE GRCh37 [hg19] v19) using STAR v2.7.10a. 