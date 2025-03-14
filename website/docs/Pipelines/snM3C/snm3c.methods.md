---
sidebar_position: 3
---
# snM3C v4.0.1 Methods

# Methods

Methylome and chromatin contact sequencing data was preprocessed for downstream analysis using the snm3C v4.0.1 pipeline (RRID:SCR_025041). Briefly, [Cutadapt software](https://cutadapt.readthedocs.io/en/stable/) was used to demultiplex paired-end sequencing reads from a single 384-well plate to cell-level FASTQ files based on a list of random primer indices, and then further used to sort, filter, and trim reads. Paired-end reads were then aligned to the human hg38 v43 reference genome using HISAT-3N. Custom python scripts from the [CEMBA GitHub repository](https://github.com/DingWB/cemba_data) were then called to separate unmapped reads, unique reads, and multi-mapped reads. The unmapped reads were saved to a FASTQ file and used for single-end alignment with HISAT-3N. Overlapping reads were removed and all resulting aligned reads merged into a single BAM. All mapped reads were deduplicated using samtools and Picard. The resulting BAM was used as input to a custom CEMBA python script for chromatin contact calling based on a 2,500 base pair threshold and as input to the [ALLCools software](https://lhqing.github.io/ALLCools/intro.html) for methylation site calling. Key summary statistics for read trimming, mapping, deduplication and chromatin contacts were then calculated and exported to a summary metrics file.

Further details regarding tools, parameters, and references used in the pipeline are available in the [YAP documentation](https://hq-1.gitbook.io/mc).
