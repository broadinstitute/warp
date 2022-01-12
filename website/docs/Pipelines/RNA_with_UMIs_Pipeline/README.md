---
sidebar_position: 1
---

# RNA with UMIs Overview

| Pipeline Version | Date Updated | Documentation Authors | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [RNAWithUMIsPipeline_v0.1.0](https://github.com/broadinstitute/warp/releases?q=RNAwithUMIs&expanded=true) | January, 2022 | [Elizabeth Kiernan](mailto:ekiernan@broadinstitute.org) & [Kaylee Mathews](mailto:kmathews@broadinstitute.org)| Please file GitHub issues in warp or contact [Kylee Degatano](mailto:kdegatano@broadinstitute.org) |

<!--- Pipeline diagram will go here --->

## Introduction to the RNA with UMIs workflow

The [RNA with UMIs Pipeline](https://github.com/broadinstitute/warp/blob/develop/pipelines/broad/rna_seq/RNAWithUMIsPipeline.wdl) is an open-source, cloud-optimized workflow for processing total RNA isolated with the Transcriptome Capture (TCap) method. TCap is a technique that hybridizes exome baits to cDNA library preparations to better facilitate RNA-sequencing in low-input or poor quality (degraded) samples. These libraries may additionally be prepared with Unique Molecular Identifiers (UMIs) which can help distinguish biological signal from noise resulting from PCR amplification. 

Overall, the workflow performs UMI correction, aligns reads to the genome, quantifies gene counts, and calculates quality metrics. The workflow produces genome- and transcriptome-aligned BAMs with indices and a merged quality metrics file. 

While this workflow was created to be used with TCap RNA-seq data, it can be used to process any bulk RNA-seq data. 

<!--- tip for methods section will go here --->

<!--- quickstart table will go here --->

## Set-up

### Installation

To download the latest release of the RNA with UMIs pipeline, see the release tags prefixed with "RNAwithUMIs" on the WARP [releases page](https://github.com/broadinstitute/warp/releases). All releases of the RNA with UMIs pipeline are documented in the [RNA with UMIs changelog](https://github.com/broadinstitute/warp/blob/develop/pipelines/broad/rna_seq/RNAWithUMIsPipeline.changelog.md). 

To search releases of this and other pipelines, use the WARP command-line tool [Wreleaser](https://github.com/broadinstitute/warp/tree/develop/wreleaser).

The RNA with UMIs pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH compliant, flexible workflow management system that supports multiple computing platforms. The workflow can also be run in [Terra](https://app.terra.bio), a cloud-based analysis platform. <!--- link to public workspace will go here --->

### Inputs

The RNA with UMIs workflow inputs are specified in JSON configuration files. Example configuration files can be found in the [test_inputs](https://github.com/broadinstitute/warp/tree/develop/pipelines/broad/rna_seq/test_inputs) folder in the WARP repository.

#### Input descriptions

The workflow takes in either a set of paired-end FASTQ files or a read group unmapped BAM that has gone through base calling. 

| Input variable name | Description | Type |
| --- | --- | --- |
| bam | Read group-specific unmapped BAM file; alternatively, paired-end FASTQ files (`r1_fastq` and `r2_fastq`) may be used. | File |
| r1_fastq | Read 1 FASTQ file; alternatively, the unmapped bam file (`bam`) may be used as input. | File |
| r2_fastq | Read 2 FASTQ file; alternatively, the unmapped bam file (`bam`) may be used as input. | File |
| read1Structure | String describing how the bases in a sequencing run should be allocated into logical reads for read 1 by fgbio's [ExtractUmisFromBam](http://fulcrumgenomics.github.io/fgbio/tools/latest/ExtractUmisFromBam.html) tool; for more information about read structures, see the [fgbio documentation](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures).  | String |
| read2Structure | String describing how the bases in a sequencing run should be allocated into logical reads for read 2 by fgbio's [ExtractUmisFromBam](http://fulcrumgenomics.github.io/fgbio/tools/latest/ExtractUmisFromBam.html) tool; for more information about read structures, see the [fgbio documentation](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures).  | String |
| starIndex | TAR file containing genome indices used for the [STAR aligner](https://github.com/alexdobin/STAR/tree/2.6.1a) | File | 
| output_basename | String used as a prefix in workflow output files | String |
| gtf | The gene annotation file (GTF) used for the [RNA-SeQC](https://github.com/getzlab/rnaseqc) tool | File | 
| platform | String used to describe the sequencing platform; only required when using FASTQ files as input. | String |
| library_name | String used to describe the library; only required when using FASTQ files as input. | String |
| platform_unit | String used to describe the platform unit; only required when using FASTQ files as input. | String |
| read_group_name | String used to describe the read group name; only required when using FASTQ files as input. | String |
| sequencing_center | String used to describe the sequencing center; only required when using FASTQ files as input; default is set to “BI”. |  String |
| ref | FASTA file used for metric collection with [Picard](https://broadinstitute.github.io/picard/) tools | File |
| refIndex | FASTA index file used for metric collection with [Picard](https://broadinstitute.github.io/picard/) tools | File |
| refDict | Dictionary file used for metric collection with [Picard](https://broadinstitute.github.io/picard/) tools | File |
| refFlat | refFlat file used for metric collection with [Picard](https://broadinstitute.github.io/picard/) tools | File |
| ribosomalIntervals | Intervals file used for RNA metric collection with [Picard](https://broadinstitute.github.io/picard/) tools | File |
| exonBedFile | Bed file used for fragment size calculations with the [RNA-SeQC](https://github.com/getzlab/rnaseqc) tool; contains non-overlapping exons  | File |

### References

<!--- Maybe reference files/inputs should have their own table in their own section? It seems a little disjointed to have them in the table above and then have the files listed here --->

The pipeline supports both hg19 and hg38 references. The reference set consists of:
1. .fasta, .fai, and .dict files
2. genome index file (STAR)
1. GTF file (RNASeQC)
1. ribosomal interval list (Picard)
1. refFlat file (Picard)

#### FASTA, index, and dictionary files
When running the workflow with the hg38 reference, we recommend using a version without HLA, ALT, and decoy contigs. These non-primary assembly contigs lead to reduced sensitivity unless the mapper is ALT-aware (e.g. bwa-mem). STAR is not alt-contig aware, so these contigs should be removed from the reference using the [modified_ref.sh script](https://github.com/broadinstitute/hydro.gen/blob/ts_rna2/Analysis/874_twist_RNA/modified_ref/modified_ref.sh).

In contrast, the hg19 reference does not have nearly as many contigs as the hg38 reference, so the workflow can be run using the standard hg19 reference stored in Broad's [public reference bucket](https://console.cloud.google.com/storage/browser/_details/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta).

#### GTF file
Genome annotation files (GTFs) contain information about genes such as the start and end coordinates of each exon, name of the gene, and the type of the transcript (e.g. protein-coding, antisense). 

The workflow uses the GENCODE v27 GTF for hg38 and v19 for hg19. We selected the v27 because it was the version used by GTEX, and the v19 because it is the latest GENCODE version available for hg19.

#### Ribosomal interval list
The workflow ribosomal interval list and the refFlat file are used by Picard metrics calculations tools. The workflow uses a custom ribosomal interval list, based on the public hg38 ribosomal interval list, which has been modified to include mitochondrial rRNA coding genes. 

#### Additional reference resources
More information about ALT contigs, HLA, decoys, alt-contig aware mapping, may be found in the following resources:

* [Which reference genome to use](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) from Heng Li's blog which argues for removing ALT contigs to improve sensitivity (for both DNA and RNA). 
* [Reference genome components](https://gatk.broadinstitute.org/hc/en-us/articles/360041155232-Reference-Genome-Components) from the GATK blog. 
* [ALT-aware mapping](https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/GPipelineAltMap_fDG.htm) as described in Illumina DRAGEN documentation. 