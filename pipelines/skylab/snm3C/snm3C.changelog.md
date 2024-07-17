# 4.0.2
2024-07-09 (Date of Last Commit)

* Updated the snM3C wdl to run on Azure; this change does not affect the snM3C pipeline

# 4.0.1
2024-06-26 (Date of Last Commit)
* Added task to untar files and output files at cell level 

# 4.0.0
2024-03-15 (Date of Last Commit)
* Reconstructed code and merged tasks to optimize pipeline and reduce cost 

# 3.0.0
2024-02-23 (Date of Last Commit)

* Updated the snM3C docker to include the latest changes to the CEMBA repository; this impacts the scientific outputs
* Added docker as a workflow-level input
* Reverted the Hisat alignments to use the --no-repeat-index parameter

# 2.0.1
2024-2-15 (Date of Last Commit)

* Updated the snM3C task memory, disk, and CPUs

# 2.0.0
2024-2-13 (Date of Last Commit)

* Merged several tasks in snM3C.wdl to reduce the cost of running this pipeline
* Removed several final outputs from snM3C.wdl 

# 1.0.1
2024-01-31 (Date of Last Commit)

* Refactored the snM3C.wdl. The outputs of the pipeline have not been affected

# 1.0.0
2023-08-01 (Date of Last Commit)

* First release of the snM3C workflow
