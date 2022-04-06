# 1.0.5
2022-03-29 (Date of Last Commit)

* Updated ingest to TDR to use transactional updates
* Added contamination outputs to the workflow

# 1.0.4
2022-03-24 (Date of Last Commit)

* Update to use references stored in Google-hosted public buckets.
* Add gsutil to PATH in picard-private docker image

# 1.0.3
2022-03-14 (Date of Last Commit)

* Update base image for picard-private docker image
* Updated to Picard version 2.26.10 and GATK version 4.2.5.0 to address log4j vulnerabilities


# 1.0.2
2022-02-18 (Date of Last Commit)

* Updated to use publicly-accessible reference and annotation files.
* Updated ribosomal intervals to include unlocalized scaffolds in the UCSC naming convention to match our reference (and renamed the file to reflect the fact that the header is not the standard GRCh38)
* Updated the STAR command line arguments, as follows:
    * Add \"--alignEndsProtrude 20 ConcordantPair\"; to rescue the case where the insert size drops below the read length and the sequencer starts to read into the adapters.
    * Removed \"--limitSjdbInsertNsj 1200000\"; the default of 1,000,000 is sufficient.
    * Removed \"--outSAMstrandField intronMotif\", defaults to \"None\"
* Updated the RNASeQC2 insert size bed file from v26 to v34
* Slightly reduced memory and disk usage on several tasks.
* Standardized memory sizing.

# 1.0.1
2022-02-01 (Date of Last Commit)

* Updated SortSam task to be able to use call caching
* Updated STAR aligner to 2.7.10a for RNAWithUMIs pipeline
* Added optional tasks to write outputs to the Terra Data Repository
* Addressed memory usage in CheckFingerprint task to allow sufficient headroom for the VM

# 1.0.0
2022-01-18 (Date of Last Commit)

* Initial Release of BroadInternalRNAWithUMIs pipeline.
* The BroadInternalRNAWithUMIs pipeline wraps the RNAWithUMIs pipeline and performs additional steps that rely on Broad specific infrastructure.
* The RNA with UMIs pipeline is an open-source, cloud-optimized workflow created for processing total RNA isolated with the Transcriptome Capture (TCap) method, but can be used to process any bullk RNA-seq data. Overall, the workflow performs UMI correction, aligns reads to the genome, quantifies gene counts, and calculates quality metrics to produces genome- and transcriptome-aligned BAMs, BAIs, and a merged quality metrics file.