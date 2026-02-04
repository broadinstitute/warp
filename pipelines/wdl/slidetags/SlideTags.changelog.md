# 1.0.7
2026-01-26 (Date of Last Commit)

* Moved inputs into new Google buckets. This change does not affect the outputs of the pipeline
* Added an optional input, billing_project, to specify a billing project when accessing files in a requester pays bucket. This change does not affect the outputs of the pipeline

# 1.0.6
2026-01-22 (Date of Last Commit)

* Added a new, defaulted input cellbender_memory_GB to Optimus; this does not affect the outputs of this pipeline
* Added a task level input, mem_size, to StarSoloFastq to expose memory settings; this does not affect the outputs of this pipeline

# 1.0.5
2025-10-24 (Date of Last Commit)

* Updated the positioning.wdl code to the latest version from the Macosko Lab repository to incorporate recent improvements

# 1.0.4
2025-10-03 (Date of Last Commit)

* Reorganized pipelines into the wdl  directory

# 1.0.3
2025-09-25 (Date of Last Commit)

* Updated the positioning.wdl code to the latest version from the Macosko Lab repository to incorporate recent improvements: added DropSift support, fixed various bugs. The Positioning script produces more metadata, saves new plotting output, and runs a more efficient algorithm

# 1.0.2
2025-07-31 (Date of Last Commit)

* Added a new optional input parameter, cellbender_hardware_memory_GB, to the CellBender tasks in the Optimus.wdl; this allows users to specify the hardware memory in GB for the CellBender tasks. The default value is set to 32 GB.

# 1.0.1
2025-06-20 (Date of Last Commit)

* Added reference genome/GTF headers to fragment file via new string inputs; this change does not affect this pipeline

# 1.0.0
2025-06-12 (Date of Last Commit)

* The first major version release for the Slide-Tags pipeline.