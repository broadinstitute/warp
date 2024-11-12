---
sidebar_position: 2
---

# ATAC Library Metrics Overview

The [ATAC pipeline](README.md) uses [SnapATAC2](https://github.com/kaizhang/SnapATAC2) to generate library-level metrics in CSV format.  


| Metric | Description |
| --- | --- |
| NHash_ID | A unique identifier used to track and reference the specific sample or dataset. |
| sequenced_reads | The total number of reads generated from the sequencing process, which includes both reads that are mapped and unmapped. |
| sequenced_read_pairs | The total number of read pairs (two reads per pair) generated from the sequencing process. This is typically half of the total sequenced reads if all reads are paired. |
| fraction_valid_barcode | The fraction of reads that contain a valid barcode, indicating the proportion of reads that are correctly assigned to a specific cell or sample. |
| fraction_Q30_bases_in_read_1 | The proportion of bases in Read 1 that have a Phred quality score of 30 or higher, indicating high-confidence base calls. |
| fraction_Q30_bases_in_read_2 | The proportion of bases in Read 2 that have a Phred quality score of 30 or higher, indicating high-confidence base calls. |
| number_of_cells | The estimated number of cells captured and sequenced in the experiment, based on the barcodes identified. |
| mean_raw_read_pairs_per_cell | The average number of raw read pairs associated with each cell, providing an indication of the sequencing depth per cell. |
| median_high-quality_fragments_per_cell | The median number of high-quality (e.g., confidently mapped) fragments associated with each cell, representing typical fragment quality across cells. |
| fraction of high-quality fragments in cells | The fraction of high-quality fragments that are associated with identified cells, indicating the proportion of good-quality data that is cell-associated. |
| fraction_of_transposition_events_in_peaks_in_cells | The fraction of transposition events within identified cells that occur within peaks, which are regions of accessible chromatin. |
| fraction_duplicates | The fraction of sequenced fragments that are duplicates, which can result from PCR amplification or other factors, indicating the redundancy in the sequencing data. |
| fraction_confidently_mapped | The fraction of sequenced fragments that are confidently mapped to the reference genome, indicating the proportion of reads that align well to the genome. |
| fraction_unmapped | The fraction of sequenced fragments that could not be mapped to the reference genome, which can indicate sequencing errors, contamination, or regions not covered by the reference. |
| fraction_nonnuclear | The fraction of sequenced fragments that are mapped to non-nuclear (e.g., mitochondrial or other organellar) DNA, providing insight into contamination or organellar activity. |
| fraction_fragment_in_nucleosome_free_region | The fraction of sequenced fragments that map to nucleosome-free regions, which are indicative of accessible chromatin. |
| fraction_fragment_flanking_single_nucleosome | The fraction of sequenced fragments that map to regions flanking single nucleosomes, indicating regions with partial chromatin accessibility. |
| tss_enrichment_score | A measure of the enrichment of transposition events at transcription start sites (TSS), indicating the accessibility of promoters across the genome. |
| fraction_of_high-quality_fragments_overlapping_TSS | The fraction of high-quality fragments that overlap transcription start sites (TSS), providing insight into promoter accessibility. |
| Number_of_peaks | The total number of peaks, or regions of accessible chromatin, identified in the dataset, representing potential regulatory elements. |
| fraction_of_genome_in_peaks | The fraction of the genome that is covered by identified peaks, indicating the extent of chromatin accessibility across the genome. |
| fraction_of_high-quality_fragments_overlapping_peaks | The fraction of high-quality fragments that overlap with identified peaks, providing an indication of the efficiency of the assay in capturing accessible regions. |
| atac_percent_target | Percent of cells recovered; value is calculated as estimated_cells/expected_cells. |


