# 1.0.17
2022-09-07 (Date of Last Commit)

* Update TDR ingest script task and docker to remove staging bucket, specify timestamp fields, and use merge ingest strategy
* Remove transcriptome bam index from output
* Updated task IngestOutputsToTDR in InternalTasks.wdl with new docker tag to accommodate changes for BroadInternalArrays pipeline. Change has no effect on this pipeline.

# 1.0.16
2022-07-15 (Date of Last Commit)

* Updated task MakeOptionalOutputBam in Utilities.wdl, this update has no effect on this pipeline

# 1.0.15
2022-06-21 (Date of Last Commit)

* Changed QC.CheckFingerprint to QC.CheckFingerprintTask to avoid a naming conflict in the update scala tests, no effect on this pipeline

# 1.0.14
2022-06-07 (Date of Last Commit)

* Fixing whitespace in the BroadInternalRNAWithUMIS.wdl, this has no functional effect on the pipeline

# 1.0.13
2022-06-03 (Date of Last Commit)

* Renamed the CompareVCFs task in VerifyCheckFingerprint.wdl to CompareVcfsAllowingQualityDifferences, this update has no effect on this pipeline

# 1.0.12
2022-06-03 (Date of Last Commit)

* Updated whitespace in BroadInternalRNAWithUMIS.wdl, this has no functional effect on the pipeline

# 1.0.11
2022-05-19 (Date of Last Commit)

* Patch security vulnerability in arrays-picard-private docker image
* Update arrays internal tasks, this update has no effect on this pipeline

# 1.0.10
2022-04-26 (Date of Last Commit)
* Remove rounding on some metrics outputs in RNAWithUMIsTasks.formatPipelineOutputs for TDR inputs
* Handle missing file inputs to TDR

# 1.0.9
2022-04-21 (Date of Last Commit)
* Update base image for picard-private docker image
* Updated to Picard version 2.26.11 and GATK version 4.2.6.1 to address log4j vulnerabilities

# 1.0.8
2022-04-20 (Date of Last Commit)

* Added memory unit to the PostprocessTranscriptomeForRSEM task in the RNAWithUMIsTasks.wdl

# 1.0.7
2022-04-12 (Date of Last Commit)

# 1.0.7
2022-04-14 (Date of Last Commit)
* Clip adapter bases pre-alignment & associated updates for TDR ingest

# 1.0.6
2022-04-04 (Date of Last Commit)

* Update task CopyFilesFromCloudToCloud in Utilities.wdl, this update has no effect on this pipeline

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