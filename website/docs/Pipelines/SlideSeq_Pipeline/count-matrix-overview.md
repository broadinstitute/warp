---
sidebar_position: 2
---

# Slide-seq Count Matrix Overview

:::warning
The Loom matrix is deprecated and the default matrix is now h5ad.
:::

The Slide-seq pipeline's default count matrix output is a h5ad file generated using [AnnData](https://anndata.readthedocs.io/en/latest/index.html).  

It contains the raw bead-by-gene counts, which vary depending on the workflow's `count_exons` parameter. By default, `count_exons` is set to `true` and the output h5ad file will contain whole-gene counts with exon counts in an additional layer. 

If the workflow is run with `count_exons` set to `false`, the output h5ad file will contain whole-gene counts. Running the workflow in this configuration will cause the h5ad matrix to have fewer columns (bead barcodes) due to the difference in STARsolo counting mode.

You can determine which type of counts are in the h5ad file by looking at the unstructured metadata (the `anndata.uns` property of the matrix) `expression_data_type` key (see [Table 1](#table-1-global-attributes) below).

The matrix also contains multiple metrics for both individual bead barcodes (the `anndata.obs` property of the matrix; [Table 2](#table-2-cell-metrics)) and individual genes (the `anndata.var` property of the matrix; [Table 3](#table-3-gene-metrics)) 

## Table 1. Global attributes

The global attributes (unstuctured metadata) in the h5ad apply to the whole file, not any specific part. 

| Attribute | Details |
| :-------- | :------ |
| `expression_data_type` | String describing if the pipeline counted whole transcript (exonic and intronic) or only exonic reads determined by the value of the `count_exons` parameter. By default, `count_exons` is `true` and `expression_data_type` is `whole_transcript`; if `count_exons` is `false` then `expression_data_type` is `exonic`. |
| `input_id` | The `input_id` provided to the pipeline as input and listed in the pipeline configuration file. This can be any string, but it's recommended for this to be consistent with any sample metadata. |
| `optimus_output_schema_version` | h5ad file spec version used during creation of the h5ad file. |
| `pipeline_version` | Version of the Slide-seq pipeline used to generate the h5ad file. |

## Table 2. Column attributes (bead barcode metrics)

The bead barcode metrics below are computed using [TagSort](https://github.com/broadinstitute/warp-tools/tree/develop/TagSort) from the [warp-tools repository](https://github.com/broadinstitute/warp-tools/tree/develop), with the exception of `input_id` which is an input to the pipeline.

| Bead Barcode Metrics | Details |
| :------------------- | :------ |
|`cell_names` | The unique identifier for each bead based on bead barcodes; identical to `CellID`. |
| `CellID` | The unique identifier for each bead based on bead barcodes; identical to `cell_names`. |
|`n_reads`| The number of reads associated with this entity. n_reads, like all metrics, are calculated from the Slide-Seq output BAM. Prior to alignment with STARsolo, reads are checked against the whitelist (1 hamming distance). These CB-corrected reads are the input to the STAR aligner. Then, the reads also get CB correction during STAR. For this reason, almost all reads in the aligned BAM have a CB tag and UB tag. Therefore, n_reads represents CB corrected reads, not all reads in the input FASTQ files. |
|`noise_reads`| Number of reads that are categorized by 10x Genomics Cell Ranger as "noise". Refers to long polymers, or reads with high numbers of N (ambiguous) nucleotides. |
|`perfect_molecule_barcodes`| The number of reads whose molecule barcodes contain no errors. |
| `reads_mapped_exonic` | The number of unique reads counted as exon; counted when BAM file's `sF` tag is assigned to `1` or `3` and the `NH:i` tag is `1`. |
| `reads_mapped_exonic_as` | The number of reads counted as exon in the antisense direction; counted when the BAM file's `sF` tag is assigned to a `2` or `4` and the `NH:i` tag is `1`. |
| `reads_mapped_intronic` | The number of reads counted as intron; counted when the BAM file's `sF` tag is assigned to a `5` and the `NH:i` tag is `1`. | 
| `reads_mapped_intronic_as` | The number of reads counted as intron in the antisense direction; counted when the BAM file's `sF` tag is assigned to a `6` and the `NH:i` tag is `1`. |
|`reads_mapped_uniquely`| The number of reads mapped to a single unambiguous location in the genome. |
|`reads_mapped_multiple`| The number of reads mapped to multiple genomic positions with equal confidence. |
| `duplicate_reads` | The number of duplicate reads. |
|`spliced_reads`| The number of reads that overlap splicing junctions. |
|`antisense_reads`| The number of reads that are mapped to the antisense strand instead of the transcribed strand. |
|`n_molecules`| Number of molecules corresponding to this entity (only reflects reads with CB and UB tags). |
|`n_fragments`| Number of fragments corresponding to this entity. |
|`fragments_with_single_read_evidence`| The number of fragments associated with this entity that are observed by only one read. |
|`molecules_with_single_read_evidence`| The number of molecules associated with this entity that are observed by only one read. |
|`perfect_cell_barcodes`| The number of reads whose bead barcodes contain no errors. |
| `reads_mapped_intergenic` | The number of reads counted as intergenic; counted when the BAM file's `sF` tag is assigned to a `7` and the `NH:i` tag is `1`. |
| `reads_unmapped` | The total number of reads that are unmapped; counted when the BAM file's `sF` tag is `0`. |
|`reads_mapped_too_many_loci`| The number of reads that were mapped to too many loci across the genome and as a consequence, are reported unmapped by the aligner. |
| `n_genes` | The number of genes detected by this bead. |
| `genes_detected_multiple_observations` | The number of genes that are observed by more than one read in this entity. |
| `molecule_barcode_fraction_bases_above_30_mean` | The average fraction of bases in molecule barcodes that receive quality scores greater than 30 across the reads of this entity. |
| `molecule_barcode_fraction_bases_above_30_variance` | The variance in the fraction of bases in molecule barcodes that receive quality scores greater than 30 across the reads of this entity.|
| `genomic_reads_fraction_bases_quality_above_30_mean` | The average fraction of bases in the genomic read that receive quality scores greater than 30 across the reads of this entity. |
| `genomic_reads_fraction_bases_quality_above_30_variance` | The variance in the fraction of bases in the genomic read that receive quality scores greater than 30 across the reads of this entity. |
| `genomic_read_quality_mean` | Average quality of base calls in the genomic reads corresponding to this entity. |
| `genomic_read_quality_variance` | Variance in quality of base calls in the genomic reads corresponding to this entity. |
|`reads_per_molecule`| The average number of reads associated with each molecule in this entity. |
|`reads_per_fragment`| The average number of reads associated with each fragment in this entity. |
|`fragments_per_molecule`| The average number of fragments associated with each molecule in this entity. |
|`cell_barcode_fraction_bases_above_30_mean`| The average fraction of base calls for the bead barcode sequences that are greater than 30, across molecules. |
|`cell_barcode_fraction_bases_above_30_variance`| The variance of the fraction of  base calls for the bead barcode sequences that are greater than 30, across molecules. |
|`n_mitochondrial_genes`| The number of mitochondrial genes detected by this bead. |
|`n_mitochondrial_molecules`| The number of molecules from mitochondrial genes detected for this bead. |
|`pct_mitochondrial_molecules`| The percentage of molecules from mitochondrial genes detected for this bead. |
| `input_id` | The `input_id` provided to the pipeline as input and listed in the pipeline configuration file. This can be any string, but it's recommended for this to be consistent with any sample metadata. |


## Table 3. Row attributes (gene metrics)

The gene metrics below are computed using [TagSort](https://github.com/broadinstitute/warp-tools/tree/develop/TagSort) from the [warp-tools repository](https://github.com/broadinstitute/warp-tools/tree/develop) except where specified.

| Gene Metrics | Details |
| ------------ | ------- |
|`gene_names` | The unique `gene_name` provided in the [GENCODE GTF](https://www.gencodegenes.org/); identical to the `Gene` attribute. |
|`ensembl_ids` | The `gene_id` provided in the [GENCODE GTF](https://www.gencodegenes.org/). |
| `Gene` | The unique `gene_name` provided in the [GENCODE GTF](https://www.gencodegenes.org/); identical to the `gene_names` attribute. |
|`n_reads`| The number of reads associated with this entity. n_reads, like all metrics, are calculated from the Slide-Seq output BAM. Prior to alignment with STARsolo, reads are checked against the whitelist (1 hamming distance). These CB-corrected reads are the input to the STAR aligner. Then, the reads also get CB correction during STAR. For this reason, almost all reads in the aligned BAM have a CB tag and UB tag. Therefore, n_reads represents CB corrected reads, not all reads in the input FASTQ files. |
|`noise_reads`| The number of reads that are categorized by 10x Genomics Cell Ranger as "noise". Refers to long polymers, or reads with high numbers of N (ambiguous) nucleotides. |
|`perfect_molecule_barcodes`| The number of reads with molecule barcodes that have no errors. |
| `reads_mapped_exonic` | The number of unique reads counted as exon; counted when BAM file's `sF` tag is assigned to `1` or `3` and the `NH:i` tag is `1`. |
| `reads_mapped_exonic_as` | The number of reads counted as exon in the antisense direction; counted when the BAM file's `sF` tag is assigned to a `2` or `4` and the `NH:i` tag is `1`. |
| `reads_mapped_intronic` | The number of reads counted as intron; counted when the BAM file's `sF` tag is assigned to a `5` and the `NH:i` tag is `1`. | 
| `reads_mapped_intronic_as` | The number of reads counted as intron in the antisense direction; counted when the BAM file's `sF` tag is assigned to a `6` and the `NH:i` tag is `1`. |
|`reads_mapped_uniquely`| The number of reads mapped to a single unambiguous location in the genome. |
|`reads_mapped_multiple`| The number of reads mapped to multiple genomic positions with equal confidence. |
| `duplicate_reads` | The number of duplicate reads. |
|`spliced_reads`| The number of reads that overlap splicing junctions. |
|`antisense_reads`| The number of reads that are mapped to the antisense strand instead of the transcribed strand. |
|`molecule_barcode_fraction_bases_above_30_mean`| The average fraction of bases in molecule barcodes that receive quality scores greater than 30 across the reads of this entity. |
|`molecule_barcode_fraction_bases_above_30_variance`| The variance in the fraction of bases in molecule barcodes that receive quality scores greater than 30 across the reads of this entity. |
|`genomic_reads_fraction_bases_quality_above_30_mean`| The average fraction of bases in the genomic read that receive quality scores greater than 30 across the reads of this entity. |
|`genomic_reads_fraction_bases_quality_above_30_variance`| The variance in the fraction of bases in the genomic read that receive quality scores greater than 30 across the reads of this entity. |
|`genomic_read_quality_mean`| Average quality of base calls in the genomic reads corresponding to this entity. |
|`genomic_read_quality_variance`| Variance in quality of base calls in the genomic reads corresponding to this entity. |
|`n_molecules`| Number of molecules corresponding to this entity (only reflects reads with CB and UB tags). |
|`n_fragments`| Number of fragments corresponding to this entity. |
|`reads_per_molecule`| The average number of reads associated with each molecule in this entity. |
|`reads_per_fragment`|The average number of reads associated with each fragment in this entity. |
|`fragments_per_molecule`| The average number of fragments associated with each molecule in this entity. |
|`fragments_with_single_read_evidence`| The number of fragments associated with this entity that are observed by only one read. |
|`molecules_with_single_read_evidence`| The number of molecules associated with this entity that are observed by only one read. |
|`number_cells_detected_multiple`| The number of bead barcodes which observe more than one read of this gene. |
|`number_cells_expressing`| The number of bead barcodes that detect this gene. |
