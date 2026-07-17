# 0.0.4
2026-07-17 (Date of Last Commit)

* * Updated tasks to use the official bcftools-vcftools and sv-imputation-rust-tools docker images
* sv-imputation-rust-tools contains the 3 rust binaries and tasks have been updated accordingly
  to use the binaries from this image instead of building them in the pipeline

# 0.0.3
2026-07-16 (Date of Last Commit)

* Optional cpu override input for Glimpse Phase

# 0.0.2
2026-07-14 (Date of Last Commit)

* remove output_prefix input from PreProcessGVCFs call

# 0.0.1
2026-07-08 (Date of Last Commit)

* Early draft of the glimpse sv imputation pipeline (putting this here to satisfy PR checks for now)
