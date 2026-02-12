# 5.1.0
2026-02-11 (Date of Last Commit)

* refactored to cleanly isolate the MitoAnnotate task and its outputs, mitochondrial deduplication is only handled within the MitoAnnotate task
* ModifyGTF and ModifyGTFMarmoset are now their own cleanly separated tasks, and the logic to determine which one to run is handled in the main workflow rather than inside BuildStarSingleNucleus
* skip_gtf_modification is now `run_modify_gtf`, a required input with no default, so you must know if your input GTF has already been modified or not
* updated the metadata.txt output to explicitly include versions and which tasks and file modifications were run in the pipeline

# 5.0.4
2026-02-11 (Date of Last Commit)

* Fixed mitofinder logic to properly handle cases where we deduplicate mitochondrial contigs from the genome fasta

# 5.0.3
2025-12-16 (Date of Last Commit)

* Remove duplicate mitochondrial contig from the genome file if mito_accession is provided

# 5.0.2
2025-11-20 (Date of Last Commit)

* Added logic to add in gene_name attribute in GTF file in SNSS2AddIntronsToGTF task

# 5.0.1
2025-11-04 (Date of Last Commit)

* Added an optional input: skip_gtf_modification (Boolean, defaulted to false)

# 5.0.0
2025-09-30 (Date of Last Commit)

* Added the MitoAnnotate task to optionally append mitochondrial sequence and annotations via MitoFinder
  This task runs early in the pipeline and updates the genome FASTA and GTF prior to building STAR and BWA indices

# 4.2.1
2025-09-26 (Date of Last Commit)

* Modified the SNSS2AddIntronsToGTF task to support Refseq GTF files in addition to GENCODE files

# 4.2.0
2025-09-10 (Date of Last Commit)

* Modified the SNSS2AddIntronsToGTF task to create a tarred STAR index alongside the intron-modified GTF file. This is a new output of the pipeline

# 4.1.0
2025-08-20 (Date of Last Commit)

* Added the SNSS2AddIntronsToGTF task to modify human and mouse GENCODE GTF files by adding introns. This is a new output of the pipeline

# 4.0.0
2025-01-17 (Date of Last Commit)

* Updated the WDL to include a new docker version 2.1.0 which has new python scripts for handling a custom marmoset GTF input
* Updated the WDL to run new marmoset scripts if the organism input is set to marmoset

# 3.1.0
2024-11-26 (Date of Last Commit)

* Added metadata.txt file as an output to the pipeline

# 3.0.0
2023-12-06 (Date of Last Commit)

* Updated BuildIndices to use bwa-mem2

# 2.2.1
2023-11-17 (Date of Last Commit)
* Updated the modify-gtf script to make it compatible with REFSEQ and GENCODE GTFs
* Added gene_versions to resulting modified GTFs
* Remove the add-introns script from the pipeline
* Removed pipe from chrom sizes task

# 2.1.2
2023-05-02 (Date of Last Commit)
* Updated the modify-gtf and add-intron scripts in the build indices docker so that PAR genes are removed from the GTF
* Added the chromosome sizes file as a pipeline output

# 2.1.1
2023-03-07 (Date of Last Commit)
* Added the BuildBWAreference task to create BWA references in addition to STAR references

# 2.1.0
2023-02-28 (Date of Last Commit)
* Added new inputs for reference genome source, build version, and annotation version
* Added checks to make sure input genome info matches the specified GTF file

# 2.0.1
2023-01-19 (Date of Last Commit)

* Added 'Disk' to task runtime sections to support running on Azure

# 2.0.0

2022-12-20 (Date of Last Commit)

* Removed all tasks unrelated to Optimus
* Updated modify_gtf.py script to support references from Refseq and Gencode
* Updated the input files to accept gtf and fasta files directly

# 1.0.1

2022-09-21 (Date of Last Commit)

* Docker image follows our guidelines
* Changed the type of biotypes from String to File so it localizes properly
* Changed the genome_fa to use the reference’s value instead of a modified_genome_fa that didn’t exist (which STAR was looking for and was then failing)

# 1.0.0

2022-02-01 (Date of Last Commit)

* Added modify_gtf.py and new docker for single nucleus smart-seq pipeline
* Update STAR to version 2.7.10a 
* Added biotypes as an input 

# 0.1.1

2021-11-15 (Date of Last Commit)

* Updated add-introns-to-gtf.py to use python3 instead of python2.

# 0.1.0

2021-05-03 (Date of Last Commit)

* Added a task to modify gtfs and fasta files and build indices for Single Nucleus Smart-seq pipeline
* Added a docker for the BuildStarSingleNucleus task

# 0.0.1

2021-04-07 (Date of Last Commit)

* Added a star dockerfile to STAR version 2.7.8a, previously it was 2.5.3a


