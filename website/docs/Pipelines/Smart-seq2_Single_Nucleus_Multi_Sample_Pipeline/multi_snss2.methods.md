---
sidebar_position: 2
---

# Smart-seq2 Single Nucleus Multi-Sample v1.3.4 Publication Methods

Below we provide an example methods section for a publication. For the complete pipeline documentation, see the [Smart-seq2 Single Nucleus Multi-Sample Overview](./README.md).

## Methods

Data preprocessing and count matrix construction for a batch (or plate) were performed using the Smart-seq2 Single Nucleus Multi-Sample v1.3.4 Pipeline (RRID:SCR_021312) as well as Picard v.2.26.10 with default tool parameters unless otherwise specified. Genomic references are publicly available in the [Broad References](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/mm10/v0/single_nucleus?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false) Google Bucket and are also listed in the [example workflow configuration](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/smartseq2_single_nucleus_multisample/mouse_example.json) in GitHub. 

For each nucleus in the batch, paired-end FASTQ files were first trimmed to remove adapters using the fastq-mcf tool with a subsampling parameter of 200,000 reads. The trimmed FASTQ files were then aligned to the GENCODE GRCm38 mouse genome using STAR v.2.7.10a. To count the number of reads per gene, but not isoforms, the quantMode parameter was set to GeneCounts. Multi-mapped reads, and optical and PCR duplicates, were removed from the resulting aligned BAM using the Picard MarkDuplicates tool with REMOVE_DUPLICATES = true. Metrics were collected on the deduplicated BAM using Picard CollectMultipleMetrics with VALIDATION_STRINGENCY =SILENT.

Intronic and exonic alignments were counted using the featureCounts v2.0.2 tool and a custom GTF that was modified from the GENCODE M23 GTF to include intronic annotations. Alignments that overlapped an annotated intron, or overlapped a single exon/intron junction, by a minimum of 3 bp were counted once as an intron. 

To count exonic alignments, a custom python script was used to create an intermediate BAM in which all alignments that cross only one intron-exon junction were removed, leaving alignments that overlap exons as well as those that overlap more than one intron-exon junction. FeatureCounts was applied to the intermediate BAM to count exonic alignments that had a 1 bp minimum overlap with an annotated exon or crossed two or more intron-exon junctions, which could be indicative of spliced RNA. 

A custom python script was then used to combine raw intronic and exonic counts as well Picard metrics into a final cell-by-gene Loom matrix.

An example of the pipeline and outputs can be found in [Terra](https://app.terra.bio/#workspaces/warp-pipelines/Smart-seq2_Single_Nucleus_Muti-Sample) and additional documentation can be found in the [Smart-seq2 Single Nucleus Multi-Sample Overview](./README.md). 