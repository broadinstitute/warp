# aou_9.0.2
2026-07-21 (Date of Last Commit)

* Added run_chrX workflow and task inputs with default false behavior
* Added optional participant_sex_tsv input and required it when run_chrX is true
* Updated the Dataproc submit command to pass participant sex TSV to chrX runs

# aou_9.0.1
2026-06-24 (Date of Last Commit)

* Added additional spark runtime parameters that are passed to the python script

* Capped the cluster prefix

# aou_9.0.0
2026-06-22 (Date of Last Commit)

* Added first version of remove_phased_samples workflow for scattering across input MatrixTable paths
* Added Dataproc launcher task pattern aligned with phasing WDL conventions in this directory
* Added basename-driven VCF/MT output naming and exposed VCF index output URLs
