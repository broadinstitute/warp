# 5.0.0
2026-11-13 (Date of Last Commit)

### Breaking changes

* Removed the CGN context-specific ALLC outputs; the unique_reads_cgn_extraction_allc_extract_array and unique_reads_cgn_extraction_tbi_extract_array outputs and the allcools extract-allc step that produced them are no longer part of the pipeline

### Additional changes

* Added a Cytozip task that converts the deduped BAMs to cytozip .cz format using czip bam_to_cz --mode full; the .cz files are emitted as new cz_tar_array and cz_count_tar_array outputs alongside the existing ALLC outputs
* Added the deduped BAMs as a new deduped_bam_tar_array output; previously they were written to a scratch directory and discarded with the VM

# 4.1.1
2025-09-19 (Date of Last Commit)

* Resolved merge conflicts and reorganize WDL pipelines into unified directory

# 4.1.0
2025-08-11 (Date of Last Commit)

* Changed read threshold max from 10M to 6M in demulitiplexing step

# 4.0.5
2025-06-16 (Date of Last Commit)

* Switched the cromwell_root_dir variable to be compatabile with Google Batch

# 4.0.4
2024-08-06 (Date of Last Commit)

* Updated the Demultiplexing task in the snm3C wdl to flag when file/cell is empty

# 4.0.3
2024-08-06 (Date of Last Commit)

* Updated the Demultiplexing task in snm3C wdl to dynamically update the batch number based on the number of fastq files present

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
