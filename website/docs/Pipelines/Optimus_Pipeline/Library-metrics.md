---
sidebar_position: 5
---

# Optimus Library-level metrics

The following table describes the library level metrics of the produced by the Optimus workflow. These are calcuated using custom python scripts available in the warp-tools repository. The Optimus workflow aligns files in shards to parallelize computationally intensive steps. This results in multiple matrix market files and shard-level library metrics. 

To produce the library-level metrics here, the [combined_mtx.py script](https://github.com/broadinstitute/warp-tools/blob/develop/3rd-party-tools/star-merge-npz/scripts/combined_mtx.py) combines all the shard-level matrix market files into one raw mtx file. Then, STARsolo is run to filter this matrix to only those barcodes that meet STARsolo's criteria of cells (using the Emptydrops_CR parameter). This matrix is then used as input during h5ad generation, and metrics are calculated from the final h5ad using the custom [add_library_tso_doublets.py](https://github.com/broadinstitute/warp-tools/tree/develop/tools/scripts) script.


| Metric | Description |
| ---| --- |
| nhash_id | The first line of of the metrics CSV echos the NHash ID if specified in the workflow run |
| number_of_reads | Total number of reads.|
| sequencing_saturation | Proportion of unique molecular identifiers (UMIs) observed relative to the total number of possible UMIs. |
| fraction_of_unique_reads_mapped_to_genome | Fraction of unique reads that map to the genome. |
| fraction_of_unique_and_multiple_reads_mapped_to_genome| Fraction of both unique and multiple reads that map to the genome. |
| fraction_of_reads_with_Q30_bases_in_rna | Fraction of reads with base quality score ≥ Q30 in RNA sequences. |
| fraction_of_reads_with_Q30_bases_in_cb_and_umi | Fraction of reads with base quality score ≥ Q30 in cell barcode (CB) and unique molecular identifier (UMI). |
| fraction_of_reads_with_valid_barcodes | Fraction of reads with valid cell barcodes. |
| reads_mapped_antisense_to_gene | Number of reads mapped antisense to gene regions.  |
| reads_mapped_confidently_exonic | Number of reads mapped confidently to exonic regions. |
| reads_mapped_confidently_to_genome | Number of reads mapped confidently to the genome. |
| reads_mapped_confidently_to_intronic_regions | Number of reads mapped confidently to intronic regions. |
| reads_mapped_confidently_to_transcriptome | Number of reads mapped confidently to the transcriptome. |
| estimated_cells | Estimated number of cells from STARsolo using the Emptydops_CR parameter. |
| umis_in_cells | Total number of unique molecular identifiers (UMIs) in cells. |
| mean_umi_per_cell | Average number of UMIs per cell. |
| median_umi_per_cell | Median number of UMIs per cell. |
| unique_reads_in_cells_mapped_to_gene | Number of unique reads in cells mapped to genes. |
| fraction_of_unique_reads_in_cells  | Fraction of unique reads in cells. |
| mean_reads_per_cell | Average number of reads per cell. |
| median_reads_per_cell | Median number of reads per cell. |
| mean_gene_per_cell | Average number of genes per cell. |
| median_gene_per_cell  | Median number of genes per cell. |
| total_genes_unique_detected | Total number of unique genes detected.  |
| percent_target | Percentage of target cells. Calculated as: estimated_number_of_cells / barcoded_cell_sample_number_of_expected_cells |
| percent_intronic_reads | Percentage of intronic reads. Calculated as: reads_mapped_confidently_to_intronic_regions / number_of_reads |
| percent_doublets | Percentage of cells flagged as doublets based on doublet scores calculated from a modified [DoubletFinder](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6853612/) algorithm. | 
| keeper_mean_reads_per_cell | Mean reads per cell for cells with >1500 genes or nuclei with >1000 genes, and doublet_score < 0.3. |
| keeper_median_genes | Median genes per cell for cells with >1500 genes or nuclei with >1000 genes, and doublet_score < 0.3>.  |
| keeper_cells | Number of cells with >1500 genes or nuclei with >1000 genes, and doublet score < 0.3.|
| percent_keeper | Percentage of keeper cells. Calculated as: keeper_cells / estimated_cells |
| percent_usable | Percentage of usable cells. Calculated as: keeper_cells / expected_cells |
| frac_tso | Fraction of reads containing TSO sequence. Calculated as the number of reads that have 20 bp or more of TSO Sequence clipped from 5' end/ total number of reads. | 