# 0.0.9
2026-08-13 (Date of Last Commit)

* Lower `PreprocessPLs` task memory to 2 GiB and increase preemptible tries to 4

# 0.0.8
2026-08-12 (Date of Last Commit)

* update multi_level_paste_pipeline_version to 0.0.4

# 0.0.7
2026-08-11 (Date of Last Commit)

* Replace per-sample array and FOFN inputs with a single `gvcf_manifest` input that is then parsed inside the workflow.
* Emit the parsed sample count so callers can use it directly without re-counting manifest rows.

# 0.0.6
2026-08-04 (Date of Last Commit)

* Remove leading `.` from `output_prefix` in `PreprocessPLsGVCF` task call

# 0.0.5
2026-07-30 (Date of Last Commit)

* remove sample_names_map_file input and associated logic/task to simplify the workflow a little

# 0.0.4
2026-07-21 (Date of Last Commit)

* Set `noAddress` to true in tasks

# 0.0.3
2026-07-17 (Date of Last Commit)

* Updated tasks to use the official sv-imputation-rust-tools docker image
* sv-imputation-rust-tools contains the 3 rust binaries and tasks have been updated accordingly
  to use the binaries from this image instead of building them in the pipeline

# 0.0.2
2026-07-14 (Date of Last Commit)

* fix some spacing in the MapSampleNames task
* remove output_prefix workflow input

# 0.0.1
2026-07-08 (Date of Last Commit)

* Early draft of the glimpse sv imputation pipeline (putting this here to satisfy PR checks for now)
