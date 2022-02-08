# 1.0.2
2022-02-01 (Date of Last Commit)

* Updated to use publicly-accessible reference and annotation files.

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