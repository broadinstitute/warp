# Whole Genome Germline Single Sample v2.2.0 Methods
The following contains a detailed methods description outlining the pipeline’s process, software, and tools that can be modified for a publication methods section.

## Detailed Methods
 
Preprocessing and variant calling was performed using the WholeGenomeGermlineSingleSample 2.2.0 pipeline using Picard 2.23.8, GATK 4.1.8, and Samtools 1.11 with default tool parameters unless otherwise specified. All reference files are available in the public [Broad References Google Bucket](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0). The pipeline follows GATK Best Practices as previously described ([Van der Auwera & O'Connor, 2020](https://www.oreilly.com/library/view/genomics-in-the/9781491975183/)) as well as the Functional Equivalence specification ([Regier et al., 2018](https://www.nature.com/articles/s41467-018-06159-4)). 

### Pre-processing and QC
Whole genome paired-end reads in unmapped BAM (uBAM) format were first scattered to perform QC and alignment in parallel. Quality metrics were calculated using Picard CollectQualityYieldMetrics. uBAMs were converted to FASTQ using Picard SamToFastq and aligned to the Hg38 reference genome using BWA mem 0.7.15 with batch size set using -K 100000000. Metadata from the uBAMs was then merged with the aligned BAMs using Picard MergeBamAlignment with the parameters --SORT_ORDER="unsorted", allowing the data to be grouped by read name for efficient downstream marking of duplicates, and --UNMAP_CONTAMINANT_READS=true, to remove cross-species contamination.

QC metrics (base distribution by cycle, insert size metrics, mean quality by cycle, and quality score distribution) were collected for the aligned, unsorted read-groups using Picard CollectMultipleMetrics. The read-group specific aligned BAMs were then aggregated and duplicate reads were flagged using Picard MarkDuplicates assuming queryname-sorted order and the parameter --OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500, which is appropriate for patterned flowcells.

The aggregated BAM file was then sorted using Picard SortSam with coordinate sort order. The fingerprints of separate read groups were verified using Picard CrosscheckFingerprints with a LOD threshold of -20. Cross-sample contamination was checked using verifyBamID2. 

The aligned BAM was then scattered for parallelization during base recalibration. A Base Quality Score Recalibration (BQSR) table was created with GATK BaseRecalibrator using the original base qualities (under the OQ Sam tag). The model was applied using ApplyBQSR with the static-quantized-quals option used according to the [Functional Equivalence specification](https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md) ([Regier et al., 2018](https://www.nature.com/articles/s41467-018-06159-4)). Recalibrated BAM files were then merged using Picard GatherSortedBamFiles. 

QC metrics were calculated for the base-recalibrated BAM using Picard CollectMultipleMetrics. Fingerprints were verified using Picard CheckFingerprint and high duplication levels and chimerism were checked using the calculated summary metrics.  

To evaluate the coverage and performance of the whole genome sequencing experiment, the BAM was assessed using the Picard tools CollectWGSMetrics and CollectRawWgsMetrics. 

The final base-recalibrated BAM was converted to CRAM using Samtools view and validated using Picard ValidateSamFile.

### Variant Calling
Prior to variant calling, the variant calling interval list was subsetted to enable parallelization. Using the GATK PrintReads tool, the WellformedReadFilter was applied to reads. Variant calling was then applied to reads that passed the filter using GATK (v3.5) HaplotypeCaller with the parameters --max_alternate_alleles 3 (sufficient for human data),  --ERC GVCF, and --read_filter OverclippedRead (to reduce false positives resulting from contamination). The resulting GVCFs were merged using Picard MergeVcfs and the final VCF file was validated using GATK ValidateVariants. Variant metrics were calculated using Picard CollectVariantCallingMetrics. 

The pipeline’s final outputs included metrics, validation reports, an aligned CRAM with index, and a GVCF containing variant calls with an accompanying index.
