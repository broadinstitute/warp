# 0.0.6
2026-05-20 (Date of Last Commit)

* Moves the VCF header reformatting step in GlimpseLigate to its own separate task that is called from main Glimpse2LowPassImputation workflow.

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
