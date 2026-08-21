# 0.0.8
2026-08-20 (Date of Last Commit)

* Updated Glimpse2SVImputationTasks dependency (sample_ids no longer required as a input)

# 0.0.7
2026-08-17 (Date of Last Commit)

* Replace `$(nproc)` usages (`--threads` args and download parallelism) with an explicit `cpu` task input on `MergeVcfs` and `ConcatVcfs`, referenced directly in the command, so concurrency always matches the task's allocated cpu

# 0.0.6
2026-08-17 (Date of Last Commit)

* make all task disk sizes dynamic and set bootDiskSizeGb to 0

# 0.0.5
2026-08-13 (Date of Last Commit)

* Remove background heartbeat monitor logging from merge task

# 0.0.4
2026-08-12 (Date of Last Commit)

* moved concat task to a shared task wdl so the batch wdl can use it

# 0.0.3
2026-07-21 (Date of Last Commit)

* Set `noAddress` to true in tasks

# 0.0.2
2026-07-17 (Date of Last Commit)

* Updated tasks to use the official bcftools-vcftools and sv-imputation-rust-tools docker images
* sv-imputation-rust-tools contains the 3 rust binaries and tasks have been updated accordingly
  to use the binaries from this image instead of building them in the pipeline

# 0.0.1
2026-07-08 (Date of Last Commit)

* Early draft of the glimpse sv imputation pipeline (putting this here to satisfy PR checks for now)
