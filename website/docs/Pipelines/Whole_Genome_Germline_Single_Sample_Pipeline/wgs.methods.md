---
sidebar_position: 2
---

# Whole Genome Germline Single Sample v3.1.6 Methods (Default workflow)

The following contains a detailed methods description outlining the pipeline’s process, software, and tools that can be modified for a publication methods section.

## Detailed methods for the default Whole Genome Germline Single Sample workflow

Preprocessing and variant calling was performed using the WholeGenomeGermlineSingleSample v3.1.6 pipeline using Picard v2.26.10, GATK v4.2.6.1, and Samtools v1.11 with default tool parameters unless otherwise specified. All reference files are available in the public [Broad References Google Bucket](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0). The pipeline follows GATK Best Practices as previously described ([Van der Auwera & O'Connor, 2020](https://www.oreilly.com/library/view/genomics-in-the/9781491975183/)) as well as the Functional Equivalence specification ([Regier et al., 2018](https://www.nature.com/articles/s41467-018-06159-4)).

### Pre-processing and quality control metrics

Whole genome paired-end reads in unmapped BAM (uBAM) format were first scattered to perform QC and alignment in parallel. Quality metrics were calculated using Picard CollectQualityYieldMetrics. uBAMs were converted to FASTQ using Picard SamToFastq and aligned to the Hg38 reference genome using BWA mem v0.7.15 with batch size set using -K 100000000. Metadata from the uBAMs was then merged with the aligned BAMs using Picard MergeBamAlignment with the parameters --SORT_ORDER="unsorted", allowing the data to be grouped by read name for efficient downstream marking of duplicates, and --UNMAP_CONTAMINANT_READS=true, to remove cross-species contamination.

QC metrics (base distribution by cycle, insert size metrics, mean quality by cycle, and quality score distribution) were collected for the aligned, unsorted read-groups using Picard CollectMultipleMetrics. The read-group specific aligned BAMs were then aggregated and duplicate reads were flagged using Picard MarkDuplicates assuming queryname-sorted order and the parameter --OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500, which is appropriate for patterned flowcells.

The aggregated BAM file was then sorted using Picard SortSam with coordinate sort order. The fingerprints of separate read groups were verified using Picard CrosscheckFingerprints ([Javed et al., 2020](https://www.nature.com/articles/s41467-020-17453-5)) with a LOD threshold of -20. Cross-sample contamination was checked using verifyBamID2.

The aligned BAM was then scattered for parallelization during base recalibration. A Base Quality Score Recalibration (BQSR) table was created with GATK BaseRecalibrator using the original base qualities (under the OQ Sam tag). The model was applied using ApplyBQSR with the static-quantized-quals option used according to the [Functional Equivalence specification](https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md) ([Regier et al., 2018](https://www.nature.com/articles/s41467-018-06159-4)). Recalibrated BAM files were then merged using Picard GatherBamFiles.

QC metrics were calculated for the base-recalibrated BAM using Picard CollectMultipleMetrics. Fingerprints were verified using Picard CheckFingerprint and high duplication levels and chimerism were checked using the calculated summary metrics.

To evaluate the coverage and performance of the whole genome sequencing experiment, the BAM was assessed using the Picard tools CollectWGSMetrics and CollectRawWgsMetrics.

The final base-recalibrated BAM was converted to CRAM using Samtools view and validated using Picard ValidateSamFile.

### Variant calling

Prior to variant calling, the variant calling interval list was split to enable parallelization. Using the GATK PrintReads tool, the WellformedReadFilter was applied to reads. Variant calling was then applied to reads that passed the filter using GATK (v3.5) HaplotypeCaller with the parameters --max_alternate_alleles 3 (sufficient for human data),  --ERC GVCF, and --read_filter OverclippedRead (to reduce false positives resulting from contamination). The resulting GVCFs were merged using Picard MergeVcfs and then reblocked using GATK ReblockGVCF with -GQB 20 -GQB 30 -GQB 40. The final reblocked GVCF file was validated using GATK ValidateVariants. Variant metrics were calculated using Picard CollectVariantCallingMetrics.

The pipeline’s final outputs included metrics, validation reports, an aligned CRAM with index, and a reblocked GVCF containing variant calls with an accompanying index.

## Detailed methods for the Functional Equivalence mode of the Whole Genome Germline Single Sample workflow

Preprocessing and variant calling was performed using the WholeGenomeGermlineSingleSample v3.1.6 pipeline using v2.26.10, GATK v4.2.6.1, and Samtools v1.11 with default tool parameters unless otherwise specified. All reference files are available in the public [Broad References Google Bucket](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0). The pipeline is functionally equivalent (as described in GATK Support: https://gatk.broadinstitute.org/hc/en-us/articles/4410456501915) to DRAGEN v3.4.12. 

### Pre-processing and quality control metrics

Whole genome paired-end reads in unmapped BAM (uBAM) format were first scattered to perform QC and alignment in parallel. Quality metrics were calculated using Picard CollectQualityYieldMetrics. uBAMs were converted to FASTQ using Picard SamToFastq and aligned to a masked version of the hg38 reference genome (reference files available at gs://gcp-public-data--broad-references/hg38/v0/dragen_reference/) using the DRAGMAP aligner. Metadata from the uBAMs was then merged with the aligned BAMs using Picard MergeBamAlignment with the parameter --UNMAP_CONTAMINANT_READS=false to maintain functional equivalence to the DRAGEN hardware.

QC metrics (base distribution by cycle, insert size metrics, mean quality by cycle, and quality score distribution) were collected for the aligned, unsorted read-groups using Picard CollectMultipleMetrics. The read-group specific aligned BAMs were then aggregated and duplicate reads were flagged using Picard MarkDuplicates assuming queryname-sorted order and the parameter --OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500, which is appropriate for patterned flowcells.

The aggregated BAM file was then sorted using Picard SortSam with coordinate sort order. The fingerprints of separate read groups were verified using Picard CrosscheckFingerprints ([Javed et al., 2020](https://www.nature.com/articles/s41467-020-17453-5)) with a LOD threshold of -20. Cross-sample contamination was checked using verifyBamID2.

QC metrics were calculated for the BAM using Picard CollectMultipleMetrics. Fingerprints were verified using Picard CheckFingerprint and high duplication levels and chimerism were checked using the calculated summary metrics.

To evaluate the coverage and performance of the whole genome sequencing experiment, the BAM was assessed using the Picard tools CollectWGSMetrics and CollectRawWgsMetrics.

The BAM was converted to CRAM using Samtools view and validated using Picard ValidateSamFile.

### Variant calling

Prior to variant calling, the DRAGEN STR model was calibrated using the CalibrateDragstrModel tool with a [public STR table](gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.str). The variant calling interval list was then split to enable parallelization. Variant calling was applied to reads using GATK4 HaplotypeCaller with the parameters --dragen-mode --disable-spanning-event-genotyping (to maintain functional equivalence to the DRAGEN hardware), and --ERC GVCF. The resulting GVCFs were merged using Picard MergeVcfs and then reblocked using GATK ReblockGVCF with -GQB 20 -GQB 30 -GQB 40. The final reblocked GVCF file was validated using GATK ValidateVariants. Variant metrics were calculated using Picard CollectVariantCallingMetrics.

The pipeline’s final outputs included metrics, validation reports, an aligned CRAM with index, and a reblocked GVCF containing variant calls with an accompanying index.

## Previous methods documents
- [WholeGenomeGermlineSingleSample_v3.0.0](https://github.com/broadinstitute/warp/blob/WholeGenomeGermlineSingleSample_v3.0.0/website/docs/Pipelines/Whole_Genome_Germline_Single_Sample_Pipeline/wgs.methods.md)
- [WholeGenomeGermlineSingleSample_v2.5.0](https://github.com/broadinstitute/warp/blob/WholeGenomeGermlineSingleSample_v2.5.0/website/docs/Pipelines/Whole_Genome_Germline_Single_Sample_Pipeline/wgs.methods.md)
- [WholeGenomeGermlineSingleSample_v2.3.7](https://github.com/broadinstitute/warp/blob/WholeGenomeGermlineSingleSample_v2.3.7/website/docs/Pipelines/Whole_Genome_Germline_Single_Sample_Pipeline/wgs.methods.md)
