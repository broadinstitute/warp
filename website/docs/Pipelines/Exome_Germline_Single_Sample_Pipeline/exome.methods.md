---
sidebar_position: 2
---

# Exome Germline Single Sample v3.0.0 Methods

The following contains a detailed methods description outlining the pipeline’s process, software, and tools that can be modified for a publication methods section.

## Detailed Methods

Preprocessing and variant calling was performed using the ExomeGermlineSingleSample 3.0.0 pipeline using Picard 2.23.8, GATK 4.2.2.0, and Samtools 1.11 with default tool parameters unless otherwise specified. All reference files are available in the public [Broad References Google Bucket](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0). The pipeline follows GATK Best Practices as previously described ([Van der Auwera & O'Connor, 2020](https://www.oreilly.com/library/view/genomics-in-the/9781491975183/)) as well as the Functional Equivalence specification ([Regier et al., 2018](https://www.nature.com/articles/s41467-018-06159-4)).

### Pre-processing and QC

Exome paired-end reads in unmapped BAM (uBAM) format were first scattered to perform QC and alignment in parallel. Quality metrics were calculated using Picard CollectQualityYieldMetrics. uBAMs were converted to FASTQ using Picard SamToFastq and aligned to the Hg38 reference genome using BWA mem 0.7.15 with batch size set using -K 100000000. Metadata from the uBAMs was then merged with the aligned BAMs using Picard MergeBamAlignment with the parameters --SORT_ORDER="unsorted", allowing the data to be grouped by read name for efficient downstream marking of duplicates, and --UNMAP_CONTAMINANT_READS=true, to remove cross-species contamination.

QC metrics (base distribution by cycle, insert size metrics, mean quality by cycle, and quality score distribution) were collected for the aligned, unsorted read-groups using Picard CollectMultipleMetrics. The read-group specific aligned BAMs were then aggregated and duplicate reads were flagged using Picard MarkDuplicates assuming queryname-sorted order and the parameter --OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500, which is appropriate for patterned flowcells.

The aggregated BAM file was then sorted using Picard SortSam with coordinate sort order. The fingerprints of separate read groups were verified using Picard CrosscheckFingerprints ([Javed et al., 2020](https://www.nature.com/articles/s41467-020-17453-5)) with a LOD threshold of -10. Cross-sample contamination was checked using verifyBamID2, using resources which had been subset to the exome target intervals.

The aligned BAM was then scattered for parallelization during base recalibration. A Base Quality Score Recalibration (BQSR) table was created with GATK BaseRecalibrator using the original base qualities (under the OQ Sam tag). The model was applied using ApplyBQSR with the static-quantized-quals option used according to the [Functional Equivalence specification](https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md) ([Regier et al., 2018](https://www.nature.com/articles/s41467-018-06159-4)). Recalibrated BAM files were then merged using Picard GatherBamFiles.

QC metrics were calculated for the base-recalibrated BAM using interval and baits lists with Picard CollectHsMetrics. Fingerprints were verified using Picard CheckFingerprint and high duplication levels and chimerism were checked using the calculated summary metrics.

The final base-recalibrated BAM was converted to CRAM using Samtools view and validated using Picard ValidateSamFile.

### Variant Calling

Prior to variant calling, the variant calling interval list was split to enable parallelization. Variant calling was performed using GATK HaplotypeCaller with the annotation parameters -G StandardAnnotation -G StandardHCAnnotation and -G AS_StandardAnnotation. Reference block genotype quality (GQ) bands were set using  -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90. The resulting GVCFs were merged using Picard MergeVcfs and reblocked using GATK ReblockGVCF with -GQB 20 -GQB 30 -GQB 40. The final reblocked GVCF file was validated using GATK ValidateVariants. Variant metrics were calculated using Picard CollectVariantCallingMetrics.

The pipeline’s final outputs included metrics, the ValidateSamFile validation reports, an aligned CRAM with index, and a reblocked GVCF containing variant calls with an accompanying index.

## Previous methods documents
- [ExomeGermlineSingleSample_v2.4.4](https://github.com/broadinstitute/warp/blob/ExomeGermlineSingleSample_v2.6.0/website/docs/Pipelines/Exome_Germline_Single_Sample_Pipeline/exome.methods.md)