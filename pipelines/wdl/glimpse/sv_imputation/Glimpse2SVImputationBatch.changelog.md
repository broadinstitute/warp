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
