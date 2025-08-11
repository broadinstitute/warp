# 1.0.19
2025-08-11 (Date of Last Commit)

* Reorganized pipelines into the wdl directory

# 1.0.18
2024-11-04 (Date of Last Commit)

* Updated to GATK version 4.6.1.0

# 1.0.17
2024-09-06 (Date of Last Commit)

* Updated to GATK version 4.6.0.0

# 1.0.16
2023-12-18 (Date of Last Commit)

* Updated to GATK version 4.5.0.0

# 1.0.15
2023-07-27 (Date of Last Commit)

* Updated UMI tools docker to address security vulnerabilities, this has no effect on this pipeline

# 1.0.14
2022-05-24 (Date of Last Commit)

* Set meta parameter `allowNestedInputs` to allow Terra users to set input parameters for heavily nested tasks.

# 1.0.13
2022-05-02 (Date of Last Commit)

* (Imported but not called by this pipeline) In MergeMetrics task, updated the list of metrics to round to integers


# 1.0.12
2023-04-07 (Date of Last Commit)

* Improvements to determinism, along with tests passing without call caching
* Allow data with very few (or 0) reads to succeed through pipeline
* Fix edge case in rnaseqc bias calculation leading to very large metrics value

# 1.0.11
2022-12-16 (Date of Last Commit)

* Updated to GATK version 4.3.0.0

# 1.0.10
2022-12-07 (Date of Last Commit)

* (Imported but not called by this pipeline) In MergeMetrics task, convert \"?\" to \"NaN\" and round some median metrics to integers.

# 1.0.9
2022-10-11 (Date of Last Commit)

* Removed task MakeOptionalOutputBam in Utilities.wdl, this update has no effect on this pipeline
* Force task rnaseqc2 to produce an empty fragment size file when rnaseqc2 does not produce this file due to insufficient data.

# 1.0.8
2022-07-29 (Date of Last Commit)

* Specify the RSEM post-processed transcriptome bam as output
* Dynamically allocate memory in Fastp task, increase fixed memory to 8gb in RNASeQC2, and increased fixed memory to 64gb in GroupByUMI
* Remove transcriptome bam index from output
* Add monitoring script to fastp and GroupByUMI tasks during soft-launch/continuous improvement
* Add maxRestries to Fastp, GroupByUMI, and RNASeQC2. Multiplier = 2 is set elsewhere.


# 1.0.7
2022-04-26 (Date of Last Commit)

* Remove rounding on some metrics outputs in RNAWithUMIsTasks.formatPipelineOutputs for TDR inputs
* Handle missing file inputs to TDR

# 1.0.6
2022-04-21 (Date of Last Commit)

* Updated to Picard version 2.26.10 and GATK version 4.2.6.1 to address log4j vulnerabilities

# 1.0.5
2022-04-20 (Date of Last Commit)

* Added memory unit to the PostprocessTranscriptomeForRSEM task in the RNAWithUMIsTasks.wdl

# 1.0.4
2022-04-08 (Date of Last Commit)

* Clip adapter bases pre-alignment

# 1.0.3
2022-04-14 (Date of Last Commit)

* Added the contamination task to the workflow.

# 1.0.2
2022-02-18 (Date of Last Commit)

* Updated the STAR command line arguments, as follows:
    * Add \"--alignEndsProtrude 20 ConcordantPair\"; to rescue the case where the insert size drops below the read length and the sequencer starts to read into the adapters.
    * Removed \"--limitSjdbInsertNsj 1200000\"; the default of 1,000,000 is sufficient.
    * Removed \"--outSAMstrandField intronMotif\", defaults to \"None\"
* Slightly reduced memory and disk usage on several tasks.
* Standardized memory sizing.

# 1.0.1
2022-02-01 (Date of Last Commit)

* Updated SortSam task to be able to use call caching.
* Updated STAR aligner to 2.7.10a for RNAWithUMIs pipeline.
* Addressed memory usage in CheckFingerprint task to allow sufficient headroom for the VM.

# 1.0.0
2022-01-14 (Date of Last Commit)

* Initial Release of RNAWithUMIs pipeline.
* The RNA with UMIs pipeline is an open-source, cloud-optimized workflow created for processing total RNA isolated with the Transcriptome Capture (TCap) method, but can be used to process any bullk RNA-seq data. Overall, the workflow performs UMI correction, aligns reads to the genome, quantifies gene counts, and calculates quality metrics to produces genome- and transcriptome-aligned BAMs, BAIs, and a merged quality metrics file.

