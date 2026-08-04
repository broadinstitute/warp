# 0.0.8
2026-08-04 (Date of Last Commit)

* Added optional `info_filter_for_inclusion` input; when > 0.0, filters both the bubble and popped
  posteriors VCFs to exclude variants with an INFO score below the given threshold

# 0.0.7
2026-08-03 (Date of Last Commit)

* Updated Glimpse2Phase task max retries and preemptible count

# 0.0.6
2026-07-31 (Date of Last Commit)

* removed unused `chunks_tsv` field from `ChunkedPanelChromosome` struct

# 0.0.5
2026-07-23 (Date of Last Commit)

* Update workflow to process array of chromosomes instead of just one
* remove unused remapping calls/task/input that were at the end of the workflow

# 0.0.4
2026-07-21 (Date of Last Commit)

* Set `noAddress` to true in tasks

# 0.0.3
2026-07-17 (Date of Last Commit)

* Updated tasks to use the official bcftools-vcftools and sv-imputation-rust-tools docker images
* sv-imputation-rust-tools contains the 3 rust binaries and tasks have been updated accordingly 
to use the binaries from this image instead of building them in the pipeline

# 0.0.2
2026-07-16 (Date of Last Commit)

* Set cpu/threads and seed for GlimpsePhase

# 0.0.1
2026-07-08 (Date of Last Commit)

* Early draft of the glimpse sv imputation pipeline (putting this here to satisfy PR checks for now)
