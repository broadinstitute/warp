---
sidebar_position: 2
---

# Slide-seq Count Matrix Overview

The Slide-seq pipeline's default count matrix output is a Loom file generated using [Loompy v.3.0.6](http://loompy.org/). 

It contains the raw bead-by-gene counts, which vary depending on the workflow's `count_exons` parameter. By default, `count_exons` is set to `true` and the output Loom will contain whole-gene counts with exon counts in an additional layer. 

If the workflow is run with `count_exons` set to `false`, the output Loom file will contain whole-gene counts. Running the workflow in this configuration will cause the Loom matrix to have fewer columns (bead barcodes) due to the difference in STARsolo counting mode.

You can determine which type of counts are in the Loom by looking at the global attribute `expression_data_type` (see [Table 1](#table-1-global-attributes) below).

The matrix also contains multiple metrics for both individual bead barcodes (the columns of the matrix; [Table 2](#table-2-column-attributes-bead-barcode-metrics)) and individual genes (the rows of the matrix; [Table 3](#table-3-row-attributes-gene-metrics)). 

## Table 1. Global attributes

The global attributes in the Loom apply to the whole file, not any specific part. 

| Attribute | Details |
| :-------- | :------ |
| `CreationDate` | Date the Loom file was created. |
| `LOOM_SPEC_VERSION` | Loom file spec version used during creation of the Loom file. |
| `expression_data_type` | String describing if the pipeline counted whole transcript (exonic and intronic) or only exonic reads determined by the value of the `count_exons` parameter. By default, `count_exons` is `true` and `expression_data_type` is `whole_transcript`; if `count_exons` is `false` then `expression_data_type` is `exonic`. |
| `input_id` | The `input_id` provided to the pipeline as input and listed in the pipeline configuration file. This can be any string, but it's recommended for this to be consistent with any sample metadata. |
| `optimus_output_schema_version` | Loom file spec version used during creation of the Loom file. |
| `pipeline_version` | Version of the Slide-seq pipeline used to generate the Loom file. |

## Table 2. Column attributes (bead barcode metrics)

The bead barcode metrics below are computed using [TagSort](https://github.com/broadinstitute/warp-tools/tree/develop/TagSort) from the [warp-tools repository](https://github.com/broadinstitute/warp-tools/tree/develop), with the exception of `input_id` which is an input to the pipeline.

| Bead Barcode Metrics | Details |
| :------------------- | :------ |
| `CellID` | The unique identifier for each bead based on bead barcodes; identical to `cell_names`. |
|`antisense_reads`| The number of reads that are mapped to the antisense strand instead of the transcribed strand. |
|`cell_barcode_fraction_bases_above_30_mean`| The average fraction of base calls for the bead barcode sequences that are greater than 30, across molecules. |
|`cell_barcode_fraction_bases_above_30_variance`| The variance of the fraction of  base calls for the bead barcode sequences that are greater than 30, across molecules. |
|`cell_names` | The unique identifier for each bead based on bead barcodes; identical to `CellID`. |
|`duplicate_reads`| The number of reads that are duplicates. |
|`fragments_per_molecule`| The average number of fragments associated with each molecule in this entity. |
|`fragments_with_single_read_evidence`| The number of fragments associated with this entity that are observed by only one read. |
| `genes_detected_multiple_observations` | The number of genes that are observed by more than one read in this entity. |
| `genomic_read_quality_mean` | Average quality of base calls in the genomic reads corresponding to this entity. |
| `genomic_read_quality_variance` | Variance in quality of base calls in the genomic reads corresponding to this entity. |
| `genomic_reads_fraction_bases_quality_above_30_mean` | The average fraction of bases in the genomic read that receive quality scores greater than 30 across the reads of this entity. |
| `genomic_reads_fraction_bases_quality_above_30_variance` | The variance in the fraction of bases in the genomic read that receive quality scores greater than 30 across the reads of this entity. |
| `input_id` | The `input_id` provided to the pipeline as input and listed in the pipeline configuration file. This can be any string, but it's recommended for this to be consistent with any sample metadata. |
| `molecule_barcode_fraction_bases_above_30_mean` | The average fraction of bases in molecule barcodes that receive quality scores greater than 30 across the reads of this entity. |
| `molecule_barcode_fraction_bases_above_30_variance` | The variance in the fraction of bases in molecule barcodes that receive quality scores greater than 30 across the reads of this entity.|
|`molecules_with_single_read_evidence`| The number of molecules associated with this entity that are observed by only one read. |
|`n_fragments`| Number of fragments corresponding to this entity. |
| `n_genes` | The number of genes detected by this bead. |
|`n_mitochondrial_genes`| The number of mitochondrial genes detected by this bead. |
|`n_mitochondrial_molecules`| The number of molecules from mitochondrial genes detected for this bead. |
|`n_molecules`| Number of molecules corresponding to this entity. |
|`n_reads`| The number of reads associated with this entity. |
|`noise_reads`| Number of reads that are categorized by 10x Genomics Cell Ranger as "noise". Refers to long polymers, or reads with high numbers of N (ambiguous) nucleotides. |
|`pct_mitochondrial_molecules`| The percentage of molecules from mitochondrial genes detected for this bead. |
|`perfect_cell_barcodes`| The number of reads whose bead barcodes contain no errors. |
|`perfect_molecule_barcodes`| The number of reads whose molecule barcodes contain no errors. |
|`reads_mapped_multiple`| The number of reads mapped to multiple genomic positions with equal confidence. |
|`reads_mapped_too_many_loci`| The number of reads that were mapped to too many loci across the genome and as a consequence, are reported unmapped by the aligner. |
|`reads_mapped_uniquely`| The number of reads mapped to a single unambiguous location in the genome. |
|`reads_per_fragment`| The average number of reads associated with each fragment in this entity. |
|`spliced_reads`| The number of reads that overlap splicing junctions. |


## Table 3. Row attributes (gene metrics)

The gene metrics below are computed using [TagSort](https://github.com/broadinstitute/warp-tools/tree/develop/TagSort) from the [warp-tools repository](https://github.com/broadinstitute/warp-tools/tree/develop) except where specified.

| Gene Metrics | Details |
| ------------ | ------- |
|`ensembl_ids` | [GENCODE GTF](https://www.gencodegenes.org/) | The gene_id listed in the GENCODE GTF. |
| `Gene` | [GENCODE GTF](https://www.gencodegenes.org/) | The unique gene_name provided in the GENCODE GTF; identical to the `gene_names` attribute. |
|`gene_names` | [GENCODE GTF](https://www.gencodegenes.org/) | The unique gene_name provided in the GENCODE GTF; identical to the `Gene` attribute. |
|`n_reads`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| The number of reads associated with this entity. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.n_reads)|
|`noise_reads`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| The number of reads that are categorized by 10x Genomics Cell Ranger as "noise". Refers to long polymers, or reads with high numbers of N (ambiguous) nucleotides. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.noise_reads)|
|`perfect_molecule_barcodes`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| The number of reads with molecule barcodes that have no errors. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.perfect_molecule_barcodes)|
|`reads_mapped_exonic`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| The number of reads for this entity that are mapped to exons. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.reads_mapped_exonic)|
|`reads_mapped_intronic`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| The number of reads for this entity that are mapped to introns. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.reads_mapped_intronic)|
|`reads_mapped_utr`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| The number of reads for this entity that are mapped to 3' untranslated regions (UTRs). [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.reads_mapped_utr)|
|`reads_mapped_uniquely`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| The number of reads mapped to a single unambiguous location in the genome. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.reads_mapped_uniquely)|
|`reads_mapped_multiple`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)|The number of reads mapped to multiple genomic positions with equal confidence. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.reads_mapped_multiple)|
|`duplicate_reads`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| The number of reads that are duplicates (see README.md for definition of a duplicate).  [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.duplicate_reads)|
|`spliced_reads`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| The number of reads that overlap splicing junctions. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.spliced_reads)|
|`antisense_reads`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| The number of reads that are mapped to the antisense strand instead of the transcribed strand. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.antisense_reads)|
|`molecule_barcode_fraction_bases_above_30_mean`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| The average fraction of bases in molecule barcodes that receive quality scores greater than 30 across the reads of this entity. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.molecule_barcode_fraction_bases_above_30_mean)|
|`molecule_barcode_fraction_bases_above_30_variance`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| The variance in the fraction of bases in molecule barcodes that receive quality scores greater than 30 across the reads of this entity. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.molecule_barcode_fraction_bases_above_30_variance)|
|`genomic_reads_fraction_bases_quality_above_30_mean`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| The average fraction of bases in the genomic read that receive quality scores greater than 30 across the reads of this entity (included for 10x Cell Ranger count comparison). [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.genomic_reads_fraction_bases_quality_above_30_mean)|
|`genomic_reads_fraction_bases_quality_above_30_variance`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| The variance in the fraction of bases in the genomic read that receive quality scores greater than 30 across the reads of this entity (included for 10x Cell Ranger count comparison). [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.genomic_reads_fraction_bases_quality_above_30_variance)|
|`genomic_read_quality_mean`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| Average quality of Illumina base calls in the genomic reads corresponding to this entity. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.genomic_read_quality_mean)|
|`genomic_read_quality_variance`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| Variance in quality of Illumina base calls in the genomic reads corresponding to this entity. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.genomic_read_quality_variance)|
|`n_molecules`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| Number of molecules corresponding to this entity. See README.md for the definition of a Molecule. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.n_molecules)|
|`n_fragments`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| Number of fragments corresponding to this entity. See README.md for the definition of a Fragment. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.n_fragments)|
|`reads_per_molecule`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| The average number of reads associated with each molecule in this entity. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.reads_per_molecule)|
|`reads_per_fragment`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| The average number of reads associated with each fragment in this entity. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.reads_per_fragment)|
|`fragments_per_molecule`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| The average number of fragments associated with each molecule in this entity. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.fragments_per_molecule)|
|`fragments_with_single_read_evidence`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| The number of fragments associated with this entity that are observed by only one read. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.fragments_with_single_read_evidence)|
|`molecules_with_single_read_evidence`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| The number of molecules associated with this entity that are observed by only one read. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.CellMetrics.molecules_with_single_read_evidence)|
|`number_cells_detected_multiple`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| The number of cells which observe more than one read of this gene. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.GeneMetrics.number_cells_detected_multiple)|
|`number_cells_expressing`|[SC Tools](https://github.com/HumanCellAtlas/sctools/tree/master/src/sctools/metrics)| The number of cells that detect this gene. [Metrics Definitions](https://sctools.readthedocs.io/en/latest/sctools.metrics.html#sctools.metrics.aggregator.GeneMetrics.number_cells_expressing)|