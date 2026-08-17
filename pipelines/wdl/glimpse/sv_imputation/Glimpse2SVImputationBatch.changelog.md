# 0.0.15
2026-08-17 (Date of Last Commit)

* Adds an optional pipeline_header_line input to the WDL

# 0.0.14
2026-08-17 (Date of Last Commit)

* make all task disk sizes dynamic and set bootDiskSizeGb to 0

# 0.0.13
2026-08-13 (Date of Last Commit)

* Remove preemptible tries for `GLIMPSE2Ligate` task
* Lower `PopAndMarginalizeCollisions` task memory to 6 GiB

# 0.0.12
2026-08-12 (Date of Last Commit)

* remove use of concat wdl and replace it with concat task from shared task wdl

# 0.0.11
2026-08-09 (Date of Last Commit)

* increase Glimpse2Phase task disk size by 20GB

# 0.0.10
2026-08-06 (Date of Last Commit)

* Update Glimpse2Ligate task memory to 18gb

# 0.0.9
2026-08-05 (Date of Last Commit)

* move re-headering commands from glimpsephase task to their own task
* add UpdateVcfSequenceDictionary to new reheadering task to keep contig header consistent across runs

# 0.0.8
2026-08-04 (Date of Last Commit)

* Update Glimpse2Phase task to stream input BCF file instead of localizing it

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
