## Announcing a new site for WARP documentation!

GDC Whole Genome Somatic Single Sample documentation has moved! Read more about the [GDC Whole Genome Somatic Single Sample Pipeline](https://broadinstitute.github.io/warp/documentation/Pipelines/Genomic_Data_Commons_Whole_Genome_Somatic/) on the new [WARP documentation site](https://broadinstitute.github.io/warp/)!

## Introduction to the GDC Whole Genome Somatic Single Sample Pipeline
The GDC Whole Genome Somatic Single Sample (abbreviated GDC here) pipeline is the alignment and preprocessing workflow for genomic data designed for the National Cancer Institute's [Genomic Data Commons](https://gdc.cancer.gov/about-gdc). 

A high-level overview of the pipeline in addition to tool parameters are available on the [GDC Documentation site](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/). 

Overall, the pipeline converts reads (CRAM or BAM) to FASTQ and (re)aligns them to the latest human reference genome (see the [GDC Reference Genome](#gdc-reference-genome) section below). Each read group is aligned separately. Read group alignments that belong to a single sample are then merged and duplicate reads are flagged for downstream variant calling. 


