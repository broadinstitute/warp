# 0.0.12
2026-06-26 (Date of Last Commit)

* use updated GLIMPSE2 build to fix cloud input streaming

# 0.0.11
2026-06-22 (Date of Last Commit)

* pre chunk each cram before passing it to mpileup

# 0.0.10
2026-06-12 (Date of Last Commit)

* remove unnecessary inputs from GlimpsePhase task to reduce cromwell metadata writing

# 0.0.9
2026-06-03 (Date of Last Commit)

* add seed input to glimpse phase task

# 0.0.8
2026-05-27 (Date of Last Commit)

* remove coverage metrics from glimpse phase task

# 0.0.7
2026-05-22 (Date of Last Commit)

* Moves the VCF header reformatting step in GlimpseLigate to its own separate task that is called from main Glimpse2LowPassImputation workflow.

# 0.0.6
2026-05-20 (Date of Last Commit)

* Updated wdl to use latest version of the Glimpse imputation image that is generated from GHA in warp-tools repo
* Updated CollectQCMetrics task to use mirror.gcr.io version of Hail image

# 0.0.5
2026-05-17 (Date of Last Commit)

* add glimpse_phase_cpu_override workflow input for mostly testing purposes

# 0.0.4
2026-05-14 (Date of Last Commit)

* pointed to correct default GLIMPSE docker image

# 0.0.3
2026-05-13 (Date of Last Commit)

* updates to support a newer version of the GLIMPSE docker image, which includes a fix for non-deterministic output during phasing when using 1 CPU.

# 0.0.2
2026-05-11 (Date of Last Commit)

* update GlimpsePhase task to use memory from ComputeShardsAndMemoryPerShard task

# 0.0.1
2026-04-30 (Date of Last Commit)

* initial implementation of low pass imputation batch wdl. this wdl will run the low pass imputation workflow on a batch of samples.
