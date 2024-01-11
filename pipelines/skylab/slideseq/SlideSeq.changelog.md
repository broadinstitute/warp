# 2.1.5
2024-01-11 (Date of Last Commit)
* Added the --soloMultiMappers flag as an optional input to the StarSoloFastq task in the StarAlign.wdl; this does affect the SlideSeq workflow

# 2.1.4
2024-01-05 (Date of Last Commit)
* Modified the STARsoloFastq task in the StarAlign.wdl so STARsolo can run different types of alignments in a single STARsolo command depending on the counting_mode; this does affect the SlideSeq workflow

# 2.1.3
2023-12-17 (Date of Last Commit)

* Updated the ATAC WDL for the Multiome BWAPairedEndAlignment and MergedBAM tasks; this does affect the SlideSeq workflow
  

# 2.1.2
2023-11-21 (Date of Last Commit)

* Added the latest warp-tools docker to tasks in the Metrics, FastqProcessing and H5adUtils wdls; this incorporates new input parameter for number of output fastq files to fastqprocess

# 2.1.1
2023-11-20 (Date of Last Commit)

* Added the latest warp-tools docker to the Metrics task; this allows use of REFSEQ references

# 2.1.0
2023-11-03 (Date of Last Commit)

* Updated the Metrics task so that Cell Metrics and Gene Metrics now calculate intronic, intronic_as, exonic, exonic_as, and intergenic metrics from unique reads only using the NH:i:1 tag in the BAM

# 2.0.1
2023-09-21 (Date of Last Commit)
* Added dynamic barcode orientation selection to the ATAC workflow FastqProcess task; this does not impact Slideseq

# 2.0.0
2023-08-22 (Date of Last Commit)

* Updated Slideseq pipeline to include STARsolo v2.7.11a; this impacts the UMI metrics CSV for unassigned genes
* Added sF tag to STARsolo aligner parameters
* Updated TagSort tool for Metrics task to calculate metrics based on the sF tag
* Modified H5adUtils task to include new metrics in the final Optimus h5ad; does not impact slide-seq
* Removed the Dropseq metrics task; this change does not impact the Slideseq workflow

# 1.0.10
2023-07-18 (Date of Last Commit)

* Added STARsolo v2.7.10b metric outputs as an optional pipeline output and an output of the STARalign and MergeSTAR tasks. This does not impact the Slideseq pipeline

# 1.0.9
2023-06-14 (Date of Last Commit)

* Updated the Metrics task WDL for adding a Dropseq task to Optimus; this has no impact on the Slide-seq workflow
* Updated the FastqProcessing.wdl for ATAC. This update has no impact on the SlideSeq workflow
* Updated STARsolo version to v2.7.10b for the StarsoloFastq task; this change does not impact this workflow
* Updated STARsolo argument for counting mode to GeneFull_Ex50pAS; this change does not impact this workflow 


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
