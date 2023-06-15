---
sidebar_position: 3
---

# Smart-seq2 Single Nucleus Multi-Sample Count Matrix Overview

The Smart-seq2 Single Nucleus Multi-Sample (Multi-snSS2) pipeline's default count matrix output is a Loom file, an HDF5 file generated using [Loompy v.3.0.6](http://loompy.org/). It contains the raw cell-by-gene intron and exon counts.

The matrix also contains multiple metrics for both individual cells (the columns of the matrix; [Table 2](#table-2-column-attributes-cell-metrics)) and individual genes (the rows of the matrix; [Table 3](#table-3-row-attributes-gene-metrics)).

Additional details for each metric are provided in the JAVA source code for Picard's [AlignmentSummaryMetrics](https://github.com/broadinstitute/picard/blob/4527d6de2f33f98f80e32e3acd32b09633529bd0/src/main/java/picard/analysis/AlignmentSummaryMetrics.java), [GcBiasSummaryMetrics](https://github.com/broadinstitute/picard/blob/4527d6de2f33f98f80e32e3acd32b09633529bd0/src/main/java/picard/analysis/GcBiasSummaryMetrics.java), and [DuplicationMetrics](https://github.com/broadinstitute/picard/blob/4527d6de2f33f98f80e32e3acd32b09633529bd0/src/main/java/picard/sam/DuplicationMetrics.java).

## Table 1. Global attributes

The global attributes in the Loom apply to the whole file, not any specific part. 

| Attribute | Details |
| :-- | :-- |
| `CreationDate` | Date the Loom file was created. |
| `LOOM_SPEC_VERSION` | Loom file spec version used during creation of the Loom file. |
| `batch_id` | The `batch_id` provided to the pipeline as input. |
| `pipeline_version` | Version of the Multi-snSS2 pipeline used to generate the Loom file. |

## Table 2. Column attributes (cell metrics)

The cell metrics below are computed using [Picard](https://broadinstitute.github.io/picard/), with the exception of `CellID`, `cell_names`, and `input_id` which are provided to the pipeline as input.

| Cell Metrics | Tool | Details |
| :----------- |:------- |:------- |
| `ACCUMULATION_LEVEL` | [CollectMultipleMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832753924763) | Level of metric accumulation; set to `ALL_READS` using the `--METRIC_ACCUMULATION_LEVEL` argument. |
| `ALIGNED_READS` | [CollectGcBiasMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832711132315) | Total number of aligned reads produced in a run. |
| `AT_DROPOUT` | [CollectGcBiasMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832711132315) | Percentage of misaligned reads with GC content below 50%. |
| `AVG_POS_3PRIME_SOFTCLIP_LENGTH.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Average length of soft-clipped bases at the 3' end of the first reads. |
| `AVG_POS_3PRIME_SOFTCLIP_LENGTH.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Average length of soft-clipped bases at the 3' end of all reads. |
| `AVG_POS_3PRIME_SOFTCLIP_LENGTH.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Average length of soft-clipped bases at the 3' end of the second reads. |
| `BAD_CYCLES.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of cycles with combined no-call and mismatch rates greater than or equal to 80% for the first reads. |
| `BAD_CYCLES.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of cycles with combined no-call and mismatch rates greater than or equal to 80% for all reads. |
| `BAD_CYCLES.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of cycles with combined no-call and mismatch rates greater than or equal to 80% for the second reads. |
| `CellID` | [warp-tools](https://github.com/broadinstitute/warp-tools) | Unique identifier for each cell provided to the pipeline as `input_ids`; identical to `cell_names` and `input_id`. |
| `ESTIMATED_LIBRARY_SIZE` | [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/13832748517275) | Estimated number of unique molecules in the library based on paired-end duplication. |
| `GC_DROPOUT` | [CollectGcBiasMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832711132315) | Percentage of misaligned reads with GC content above 50%. |
| `GC_NC_0_19` | [CollectGcBiasMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832711132315) | Normalized coverage over reads with GC content from 0 - 19%. |
| `GC_NC_20_39` | [CollectGcBiasMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832711132315) | Normalized coverage over reads with GC content from 20 - 39%. |
| `GC_NC_40_59` | [CollectGcBiasMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832711132315) | Normalized coverage over reads with GC content from 40 - 59%. |
| `GC_NC_60_79` | [CollectGcBiasMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832711132315) | Normalized coverage over reads with GC content from 60 - 79%. |
| `GC_NC_80_100` | [CollectGcBiasMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832711132315) | Normalized coverage over reads with GC content from 80 - 100%. |
| `MAD_READ_LENGTH.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Median absolute deviation of the lengths of forward reads. |
| `MAD_READ_LENGTH.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Median absolute deviation of the lengths of all reads. |
| `MAD_READ_LENGTH.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Median absolute deviation of the lengths of reverse reads. |
| `MAX_READ_LENGTH.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Maximum length of forward reads. |
| `MAX_READ_LENGTH.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Maximum length of all reads. |
| `MAX_READ_LENGTH.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Maximum length of reverse reads. |
| `MEAN_READ_LENGTH.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Mean length of forward reads. |
| `MEAN_READ_LENGTH.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Mean length of all reads. |
| `MEAN_READ_LENGTH.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Mean length of reverse reads. |
| `MEDIAN_READ_LENGTH.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Median length of forward reads. |
| `MEDIAN_READ_LENGTH.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Median length of all reads. |
| `MEDIAN_READ_LENGTH.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Median length of reverse reads. |
| `MIN_READ_LENGTH.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Minimum length of forward reads. |
| `MIN_READ_LENGTH.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Minimum length of all reads. |
| `MIN_READ_LENGTH.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Minimum length of reverse reads. |
| `PCT_ADAPTER.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of pass-filter forward reads that are unaligned or aligned with a mapping quality of 0 and match to a known adapter sequence from the start of the read. |
| `PCT_ADAPTER.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of all pass-filter reads that are unaligned or aligned with a mapping quality of 0 and match to a known adapter sequence from the start of the read. |
| `PCT_ADAPTER.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of pass-filter reverse reads that are unaligned or aligned with a mapping quality of 0 and match to a known adapter sequence from the start of the read. |
| `PCT_CHIMERAS.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of forward reads where the insert is larger than 100 kb or the ends of the pair map to different chromosomes. |
| `PCT_CHIMERAS.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of all reads where the insert is larger than 100 kb or the ends of the pair map to different chromosomes. |
| `PCT_CHIMERAS.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of reverse reads where the insert is larger than 100 kb or the ends of the pair map to different chromosomes. |
| `PCT_HARDCLIP.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of pass-filter bases that are hard-clipped from aligned, forward reads. |
| `PCT_HARDCLIP.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of pass-filter bases that are hard-clipped from all aligned reads. |
| `PCT_HARDCLIP.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of pass-filter bases that are hard-clipped from aligned, reverse reads. |
| `PCT_PF_READS.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of forward reads that pass vendor check (pass-filter). |
| `PCT_PF_READS.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of reads that pass vendor check (pass-filter). |
| `PCT_PF_READS.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of reverse reads that pass vendor check (pass-filter). |
| `PCT_PF_READS_ALIGNED.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of pass-filter forward reads that are aligned. |
| `PCT_PF_READS_ALIGNED.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of all pass-filter reads that are aligned. |
| `PCT_PF_READS_ALIGNED.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of pass-filter reverse reads that are aligned. |
| `PCT_PF_READS_IMPROPER_PAIRS.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of forward reads not properly aligned in pairs. |
| `PCT_PF_READS_IMPROPER_PAIRS.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of reads not properly aligned in pairs. |
| `PCT_PF_READS_IMPROPER_PAIRS.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of reverse reads not properly aligned in pairs. |
| `PCT_READS_ALIGNED_IN_PAIRS.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of forward reads properly aligned in pairs. |
| `PCT_READS_ALIGNED_IN_PAIRS.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of reads properly aligned in pairs. |
| `PCT_READS_ALIGNED_IN_PAIRS.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of reverse reads properly aligned in pairs. |
| `PCT_SOFTCLIP.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of pass-filter bases that are soft-clipped from aligned forward reads. |
| `PCT_SOFTCLIP.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of pass-filter bases that are soft-clipped from all aligned reads. |
| `PCT_SOFTCLIP.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of pass-filter bases that are soft-clipped from aligned reverse reads. |
| `PERCENT_DUPLICATION` |[MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/13832748517275) | Fraction of mapped sequence marked as duplicate. |
| `PF_ALIGNED_BASES.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Total number of aligned bases in pass-filter forward reads. |
| `PF_ALIGNED_BASES.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Total number of aligned bases in all pass-filter reads. |
| `PF_ALIGNED_BASES.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Total number of aligned bases in pass-filter reverse reads. |
| `PF_HQ_ALIGNED_BASES.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of bases aligned to the reference sequence in forward reads with high mapping quality. |
| `PF_HQ_ALIGNED_BASES.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of bases aligned to the reference sequence in all reads with high mapping quality. |
| `PF_HQ_ALIGNED_BASES.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of bases aligned to the reference sequence in reverse reads with high mapping quality. |
| `PF_HQ_ALIGNED_Q20_BASES.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Subset of `PF_HQ_ALIGNED_BASES.FIRST_OF_PAIR` with a base call quality of at least 20. |
| `PF_HQ_ALIGNED_Q20_BASES.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Subset of `PF_HQ_ALIGNED_BASES.PAIR` with a base call quality of at least 20. |
| `PF_HQ_ALIGNED_Q20_BASES.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Subset of `PF_HQ_ALIGNED_BASES.SECOND_OF_PAIR` with a base call quality of at least 20. |
| `PF_HQ_ALIGNED_READS.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of pass-filter forward reads aligned with a mapping quality of at least 20. |
| `PF_HQ_ALIGNED_READS.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of all pass-filter reads aligned with a mapping quality of at least 20. |
| `PF_HQ_ALIGNED_READS.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of pass-filter reverse reads aligned with a mapping quality of at least 20. |
| `PF_HQ_ERROR_RATE.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of bases in pass-filter, high-quality forward reads that do not match the reference. |
| `PF_HQ_ERROR_RATE.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of bases in all pass-filter, high-quality reads that do not match the reference. |
| `PF_HQ_ERROR_RATE.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Fraction of bases in pass-filter, high-quality reverse reads that do not match the reference. |
| `PF_HQ_MEDIAN_MISMATCHES.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Median number of mismatches in high-quality forward reads. |
| `PF_HQ_MEDIAN_MISMATCHES.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Median number of mismatches in all high-quality reads. |
| `PF_HQ_MEDIAN_MISMATCHES.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Median number of mismatches in high-quality reverse reads. |
| `PF_INDEL_RATE.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of insertion and deletion events per 100 aligned bases in forward reads. |
| `PF_INDEL_RATE.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of insertion and deletion events per 100 aligned bases in all reads. |
| `PF_INDEL_RATE.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of insertion and deletion events per 100 aligned bases in reverse reads. |
| `PF_MISMATCH_RATE.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Rate of base mismatching for all aligned bases in forward reads. |
| `PF_MISMATCH_RATE.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Rate of base mismatching for all aligned bases in all reads. |
| `PF_MISMATCH_RATE.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Rate of base mismatching for all aligned bases in reverse reads. |
| `PF_NOISE_READS.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of pass-filter forward reads marked as noise. |
| `PF_NOISE_READS.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of all pass-filter reads marked as noise. |
| `PF_NOISE_READS.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of pass-filter reverse reads marked as noise. |
| `PF_READS.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of forward reads that pass vendor check (pass-filter). |
| `PF_READS.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of reads that pass vendor check (pass-filter). |
| `PF_READS.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of reverse reads that pass vendor check (pass-filter). |
| `PF_READS_ALIGNED.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of pass-filter forward reads that are aligned. |
| `PF_READS_ALIGNED.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of pass-filter reads that are aligned. |
| `PF_READS_ALIGNED.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of pass-filter reverse reads that are aligned. |
| `PF_READS_IMPROPER_PAIRS.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of forward reads not properly aligned in pairs. |
| `PF_READS_IMPROPER_PAIRS.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of reads not properly aligned in pairs. |
| `PF_READS_IMPROPER_PAIRS.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of reverse reads not properly aligned in pairs. |
| `READS_ALIGNED_IN_PAIRS.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of forward reads properly aligned in pairs. |
| `READS_ALIGNED_IN_PAIRS.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of reads properly aligned in pairs. |
| `READS_ALIGNED_IN_PAIRS.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of reverse reads properly aligned in pairs. |
| `READS_USED` | [CollectGcBiasMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832711132315) | String describing whether duplicates are included in metrics produced by [CollectGcBiasMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832711132315); the pipeline removes duplicates before metrics are calculated. |
| `READ_PAIRS_EXAMINED` | [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/13832748517275) | Number of mapped read pairs examined by [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/13832748517275). |
| `READ_PAIR_DUPLICATES` | [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/13832748517275) | Number of read pairs marked as duplicates. |
| `READ_PAIR_OPTICAL_DUPLICATES` | [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/13832748517275) | Number of read pairs duplicates caused by optical duplication. |
| `SD_READ_LENGTH.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Standard deviation of forward read lengths. |
| `SD_READ_LENGTH.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Standard deviation of read lengths. |
| `SD_READ_LENGTH.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Standard deviation of reverse read lengths. |
| `SECONDARY_OR_SUPPLEMENTARY_RDS` | [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/13832748517275) | Number of secondary or supplemetary reads. |
| `STRAND_BALANCE.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of pass-filter forward reads aligned divided by the total number of pass-filter reads aligned. |
| `STRAND_BALANCE.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Average strand balance of forward and reverse reads. |
| `STRAND_BALANCE.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Number of pass-filter reverse reads aligned divided by the total number of pass-filter reads aligned. |
| `TOTAL_CLUSTERS` | [CollectGcBiasMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832711132315) | Total number of reads after filtering used in GC bias calculation. |
| `TOTAL_READS.FIRST_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Total number of forward reads. |
| `TOTAL_READS.PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Total number of reads. |
| `TOTAL_READS.SECOND_OF_PAIR` | [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832683660955) | Total number of reverse reads. |
| `UNMAPPED_READS` | [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/13832748517275) | Number of unmapped reads examined by [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/13832748517275). |
| `UNPAIRED_READS_EXAMINED` | [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/13832748517275) | Number of mapped reads without a mapped mate pair examined by [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/13832748517275). |
| `UNPAIRED_READ_DUPLICATES` | [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/13832748517275) | Number of fragments marked as duplicates. |
| `WINDOW_SIZE` | [CollectGcBiasMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832711132315) | Genomic window size used in GC bias calculation. |
| `cell_names` | [warp-tools](https://github.com/broadinstitute/warp-tools) | Unique identifier for each cell provided to the pipeline as `input_ids`; identical to `Cell_ID` and `input_id`. |
| `input_id` | [warp-tools](https://github.com/broadinstitute/warp-tools) | Unique identifier for each cell provided to the pipeline as `input_ids`; identical to `Cell_ID` and `cell_names`. |

## Table 3. Row attributes (gene metrics)

| Gene Metrics | Tool | Details |
| ------------ | ---- | ------- |
| `Gene` | [GENCODE GTF](https://www.gencodegenes.org/) | The unique `gene_ids` provided in the [GENCODE GTF](https://www.gencodegenes.org/); identical to the `ensembl_ids` attribute. |
| `ensembl_ids` | [GENCODE GTF](https://www.gencodegenes.org/) | The unique `gene_ids` provided in the [GENCODE GTF](https://www.gencodegenes.org/); identical to the `Gene` attribute. |
| `exon_lengths` | [warp-tools](https://github.com/broadinstitute/warp-tools) | The length in base pairs of the exons corresponding to this entity. |
| `gene_names` | [GENCODE GTF](https://www.gencodegenes.org/) | The unique `gene_name` provided in the [GENCODE GTF](https://www.gencodegenes.org/). |
| `intron_lengths` | [warp-tools](https://github.com/broadinstitute/warp-tools) | The length in base pairs of the introns corresponding to this entity. |