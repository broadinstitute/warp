# 1.0.2
2022-09-07 (Date of Last Commit)

* Increased disk space in the MakeOptionalOutputBam task in Utilities.wdl
* Updated task IngestOutputsToTDR in InternalTasks.wdl with new docker tag to accommodate changes for BroadInternalArrays pipeline. Change has no effect on this pipeline.

# 1.0.1
2022-06-21 (Date of Last Commit)

* Changed QC.CheckFingerprint to QC.CheckFingerprintTask to avoid a naming conflict in the update scala tests, no effect on this pipeline

# 1.0.0
2022-05-05 (Date of Last Commit)

* Initial Release of UltimaGenomicsWholeGenomeGermline pipeline.
* The UltimaGenomicsWholeGenomeGermline pipeline is an open-source, cloud-optimized workflow created for processing Ultima Genomics Whole Genome Sequenced Germline samples. Overall, the workflow aligns reads to the genome, marks duplicates, calls variants, and calculates quality metrics to produce a CRAM, CRAI, GVCF, filtered VCF, and quality metrics.