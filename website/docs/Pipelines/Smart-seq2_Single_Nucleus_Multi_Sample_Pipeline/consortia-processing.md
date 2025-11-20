---
sidebar_position: 4
---

# Consortia Data Processing

## Brain Initiative Cell Census Network Processing
The Smart-seq2 Single Nucleus Multi-Sample (Multi-snSS2) pipeline supports data processing for the [BRAIN Initiative Cell Census Network (BICCN)](https://biccn.org/). An overview of the BICCN pipeline resources is available on the BICCN's [Pipelines page](https://biccn.org/tools/biccn-pipelines).

### Multi-snSS2 reference files for BICCN data processing
The BICCN 2.0 Whole Mouse Brain Working Group uses the Ensembl GRCm38 reference for alignment and a modified GTF for gene annotation (see table below). All Multi-snSS2 pipeline reference inputs were created with the [BuildIndices workflow](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/build_indices).

 BICCN processes single-nucleus data, which is enriched in pre-mRNAs containing introns. To account for this, the Multi-snSS2 workflow counts reads that map to both exonic and intronic regions (any part of a contig that is not exonic nor intergenic). The BuildIndices workflow uses the `BuildStarSingleNucleus` task to add intron annotations to the GTF with a custom [python script](https://github.com/broadinstitute/warp-tools/blob/develop/3rd-party-tools/build-indices/add-introns-to-gtf.py). These annotations enable intron counting with the [featureCounts](http://subread.sourceforge.net/) software. 

 The custom GTF contains all annotations for any `gene_id` that has at least one transcript. This reduces the number of genes in the GTF to \~32,000. 

All reference files are available in a public Google bucket (see table below) and are accompanied by a README that details reference provenance (gs://gcp-public-data--broad-references/mm10/v0/README_mm10_singlecell_gencode.txt). 

| Multi-snSS2 reference input name | Google bucket URI | Reference source | Description |
| --- | --- | --- | --- |
| `annotations_gtf` | gs://gcp-public-data--broad-references/mm10/v0/single_nucleus/modified_gencode.vM23.primary_assembly.annotation.gtf | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gzf | Modified GENCODE GTF including intron annotations that can be used for intron counting with featureCounts. |
| `genome_ref_fasta` | gs://gcp-public-data--broad-references/mm10/v0/single_nucleus/modified_mm10.primary_assembly.genome.fa | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.p6.genome.fa.gz | FASTA file used to create the STAR reference files. |
| `tar_star_reference` | gs://gcp-public-data--broad-references/mm10/v0/single_nucleus/star/modified_star_2.7.9a_primary_gencode_mouse_vM23.tar | NA — built with the BuildIndices workflow. | Reference files used for alignment with STAR. |
| `adapter_list` | gs://broad-gotc-test-storage/MultiSampleSmartSeq2SingleNucleus/adapters/Illumina_adapters_list.fa | See Illumina's overview on [adapter sequences](https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html). | List of adapter sequences used for trimming. |






 









