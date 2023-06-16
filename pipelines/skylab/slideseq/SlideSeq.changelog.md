# 1.0.9
2023-06-14 (Date of Last Commit)

* Updated the Metrics task WDL for adding a Dropseq task to Optimus; this has no impact on the Slide-seq workflow
* Updated the FastqProcessing.wdl for ATAC. This update has no impact on the SlideSeq workflow 


# 1.0.8
2023-05-31 (Date of Last Commit)

* Updated the FastqProcessing.wdl. This update has no impact on the SlideSeq workflow

# 1.0.7
2023-05-11 (Date of Last Commit)

* Updated Docker image for FastqProcessing task to latest warp-tools

# 1.0.6
2023-05-04 (Date of Last Commit)

* Updated the warp-tools docker for Optimus fastqprocess; change does not impact this workflow

# 1.0.5
2023-04-23 (Date of Last Commit)

* Updated the STARalign task; does not affect this workflow

# 1.0.3
2023-04-19

* Updated warp-tools docker which included a fix for a small bug in create_snrna_optimus.py that was causing the script not to run

# 1.0.3
2023-03-27 (Date of Last Commit)
* Removed the following columns from the gene metrics csv and the Loom as the counts were empty/incorrect: reads_unmapped, reads_mapped_exonic, reads_mapped_intronic, reads_mapped_utr, reads_mapped_intergenic, duplicate_reads. We also removed duplicate_reads from the cell metrics csv and the Loom

# 1.0.2
2023-03-15 (Date of Last Commit)

* Removed the following columns from the cell metrics csv and the Loom as the counts were empty/incorrect: reads_unmapped, reads_mapped_exonic, reads_mapped_intronic, reads_mapped_utr, reads_mapped_intergenic
* Updated warp-tools docker. The latest loom building script updates the optimus_output_schema_version from 1.0.0 to 1.0.1 to capture the metrics changes listed above 


# 1.0.1
2023-02-28 (Date of Last Commit)

* Added a new task to the workflow that reads the tar_star_reference file to obtain the genomic reference source, build version, and annotation version and outputs the information as a txt file

# 1.0.0
2023-02-13 (Date of Last Commit)

* The first major version release for the SlideSeq pipeline
