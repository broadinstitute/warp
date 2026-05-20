# 0.0.16
2026-05-20 (Date of Last Commit)

* Adds an optional pipeline_header_line input to the WDL
* Adds a new task called UpdateHeader that is called to fo VCF header reformatting

# 0.0.15
2026-05-20 (Date of Last Commit)

* Updated `batch_pipeline_version` to 0.0.6 to reflect update to Glimpse and Hail docker images being used in the Batch wdl

# 0.0.14
2026-05-20 (Date of Last Commit)

* Updated `quota_consumed_version` to 0.0.3 and `input_qc_version` to 1.0.2 to reflect interface updates in auxiliary workflows

# 0.0.13
2026-05-19 (Date of Last Commit)

* Updated `ExtractAnnotations`, `RecomputeAndAnnotate`, `MergeSampleChunksVcfsWithPaste` tasks to require fewer CPUs and/or less RAM

# 0.0.12
2026-05-17 (Date of Last Commit)

* add glimpse_phase_cpu_override workflow input for mostly testing purposes

# 0.0.11
2026-05-14 (Date of Last Commit)

* Improved memory performance of RecomputeAndAnnotate to better handle large amounts of samples.
* Use correct MergeCoverageMetrics task for recombining batched outputs.

# 0.0.10
2026-05-13 (Date of Last Commit)

* Update GLIMPSE docker image to latest, which includes a fix for non-deterministic output during phasing when using 1 CPU

# 0.0.9
2026-05-11 (Date of Last Commit)

* add batch wdl pipeline version

# 0.0.8
2026-05-08 (Date of Last Commit)

* Add support for optional cram_manifest input

# 0.0.7
2026-04-30 (Date of Last Commit)

* Add batched version of the wdl, which runs n samples at a time as subworkflows and recombines their outputs.

# 0.0.6
2026-04-15 (Date of Last Commit)

* Add Input QC wdl version

# 0.0.5
2026-04-03 (Date of Last Commit)

* add quota consumed wdl version

# 0.0.4
2026-04-01 (Date of Last Commit)

* split out imputed hom ref sites to their own sites only vcf file output

# 0.0.3
2026-03-25 (Date of Last Commit)

* Reorganize wdl to be able to run on contigs more easily.  Now the workflow is fully driven by the `contigs` input
* The wdl now expects the reference related files to all live under the same cloud base path

# 0.0.2
2026-03-19 (Date of Last Commit)

* Optimizations to pipeline to run faster and more cheaply using the aou + anvil reference panel

# 0.0.1
2026-02-26 (Date of Last Commit)

* Early draft of the glimpse lowpass imputation pipeline (putting this here to satisfy PR checks for now)
