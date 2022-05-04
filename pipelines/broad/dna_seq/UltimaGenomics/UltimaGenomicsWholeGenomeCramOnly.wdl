version 1.0

import "../../../../tasks/broad/UltimaGenomicsWholeGenomeGermlineTasks.wdl" as Tasks
import "../../../../tasks/broad/Utilities.wdl" as Utilities
import "../../../../tasks/broad/GermlineVariantDiscovery.wdl" as VariantDiscoverTasks
import "../../../../tasks/broad/UltimaGenomicsWholeGenomeGermlineAlignmentMarkDuplicates.wdl" as UltimaGenomicsWholeGenomeGermlineAlignmentMarkDuplicates
import "../../../../tasks/broad/InternalTasks.wdl" as InternalTasks
import "../../../../tasks/broad/Qc.wdl" as QC
import "../../../../tasks/broad/UltimaGenomicsWholeGenomeGermlineQC.wdl" as UltimaGenomicsWholeGenomeGermlineQC
import "../../../../structs/dna_seq/UltimaGenomicsWholeGenomeGermlineStructs.wdl" as Structs

# CHANGELOG
#  1.1.1     get multiple input cram
#  1.1.2     optional gvcf, copy report outputs to Nexus
#  2.1.0     gatk 5.0 (not for real, I forgot to update the docker)
#  2.1.1     gatk 5.0 for real
#  2.1.2     changed mark duplicates disk logic
#  2.1.3     changed CheckContamination docker and updated args (--adjust-MQ 0), updated default HaplotypeCaller flags to create bam files
#  2.1.4     fixed logic in H1->H2 retry mechanism, new gatk with bugfix
#  2.1.5     new gatk/picard docker for MarkDuplicatesSpark and CrosscheckFingerprints
#  2.1.5.1   light output naming refactoring
#  2.1.5.2   MarkDuplicatesSpark disk fix
#  2.1.5.3   reduced cpu number for MarkDuplicatesSpark
#  2.1.6     minor restructuring in preparation of gatk 0.5.1
#  2.1.7     replaced ultima markduplicates docker with broad one
#  2.2.0     gatk 0.5.1, blocklist in filtering
#  2.2.1     gatk 0.5.1-2 in markduplicates without localization
#  2.2.2     added CollectIntervalCoverageStats, removed maxretries from tasks without default preemptible tries
#  2.3.1     supports V4 of the basecalling pipeline, aggregates sequencing metrics
#  2.3.2     BIOIN-73 Flow for converted libraries
#            BIOIN-80 Contig filtering replaced with allele filtering in the GATK (much lower cost for high coverage runs)
#            BIOIN-72 Expanded variant calling report
#            BIOIN-48 MarkDuplicates for single end reads
#            BIOIN-69 SV calls
#  2.3.2.1   increased memory in CreateReport
#  2.4.1     [BIOIN-91, BIOIN-114] Improved report
#            BIOIN-130 Added new ground truth datasets
#            BIOIN-109 Variant filtering now outputs correct false positive rate
#            BIOIN-49 Flow based mark duplicates turned on by default with optimized parameters
#            BIOIN-112 SOR-based allele filtering in GATK
#  2.4.2     Variant calling report without ground truth [BIOIN-103]
#            GATK updated to 4.2.0 [BIOIN-139]
#            T/N support in Mutect [BIOIN-115]
#            GATK is about 40% faster [BIOIN-138]
#            Variant report outputs added to the aggregated metrics [BIOIN-23]
#            inputs/outputs are ingested into the database
#            Funcotator reports [BIOIN-84]
# 2.5.1      Major changes
#              [BIOIN-74] Contamination is calculated from realigned reads and flow based model rather than PairHMM
#              [BIOIN-97] FeatureMapper integrated into jukebox_vc pipeline
#              [BIOIN-107, BIOIN-156] Optimization of the coverage collection code
#              [BIOIN-142] Filtering alleles without ground truth
#              [BIOIN-151] Expanded test set
#            Minor changes
#              SB and AS_HmerLength are now reported in the VCF
#              Another error code to SW crash and optional SV calls
#              Barcode and read group ID are extracted and sent to metadata
#              Changed markDuplicates parameter (expected to decrease significantly duplicate marking)
#              Changed GATK parameters (to use adaptive pruning and dynamic read disqualification)
# 2.6.1      BIOIN-168 Compatible to BARC 4.2
# 2.6.2      Major changes
#              BIOIN-142 Improved allele filtering and performance on exome
#              BIOIN-146 UA is an optional aligner
#              BIOIN-162 SNP Error rate is calculated for all runs
#            Minor changes
#              Various bug fixes and stability improvements
#               * Fixes [BIOIN-185]
#               * Fixes [BIOIN-195]
#               * Fixes [BIOIN-198]
#               * Fixes [BIOIN-202]
#               * Remedy for [BIOIN-181]
#  2.6.3     Disabled featureMap by default
#            Added exome weights
#            Fixed output parameter
#  2.7.1     [BIOIN-153] GATK outputs isolated weak variants
#            GATK rebased to version 4.2.2.0
#            [BIOIN-224] Additional metrics: RLQ30, RLQ25, MEDIAN_READ_LENGTH output by picard for the FAT table
#            [BIOIN-127] Coverage statistics outptus bigwig files, coverage annotation in comparison scripts much faste
#            Stability improvements and bug fix in the somatic pipeline
#            [BIOIN-192] Metrics from no ground truth report added to the database
#            [BIOIN-169] Removed task to remove symbolic allele
#            Minor fixes:
#             [BIOIN-183] Better output of the report metadata
#             [BIOIN-152] Improved blocklist for somatic pipeline
#  2.7.2     Allow empty adapter in converted libraries
#  3.1.1     [BIOIN-189] contamination is calcuated on the FlowBased model
#            [BIOIN-257] Apply filtering model on the gVCF
#            [BIOIN-228] VCF contain flow-based annotations
#            [BIOIN-269] Improved filteirng with no overfitting and new training sets
#            Other minor fixes
#  3.2.1     [BIOIN-288] filter reads by rsq
#            FlowBasedHMM used by default as the error model in haplotypeCaller
#            Picard overrides in metrics collection
#            Contamination uses cycle-skip files
#            Annotations in the GATK from the allele filteirng
#  3.2.2     Fixed an error when the pipeline was running without the variant calling
#            [BIOIN-260] Rewired the flow of the VCFs in the pipeline so that the final output contains interval annotation and funcotator receives the correct vcf
#  3.2.3     Removed the rsq filtering
#            Switched error model back to FlowBased

workflow UltimaGenomicsWholeGenomeCramOnly {
  input {

    ContaminationSites contamination_sites
    AlignmentReferences alignment_references
    VcfPostProcessing vcf_post_processing

    # Sample Information
    Array[File]? input_cram_list
    Array[File]? input_bam_list
    String base_file_name

    Float rsq_threshold = 1.0
    Int reads_per_split = 20000000

    Boolean save_bam_file = false
  }

  meta {
    allowNestedInputs: true
  }

  parameter_meta {
    contamination_sites: "Struct containing files for contamination estimation"
    alignment_references: "Struct containing reference files for alignment with BWA mem"
    input_cram_list: "Array of CRAM files to be used as workflow input. Must be specified if `input_bam_list` is not provided"
    input_bam_list: "Array of unmapped BAM files to be used as workflow input. Must be specified if `input_cram_list` is not provided"
    base_file_name: "Base name for each of the output files."
    rsq_threshold: "Threshold for a read quality metric that is produced by the sequencing platform"
    reads_per_split: "Number of reads by which to split the CRAM prior to alignment"
    save_bam_file: "If true, then save intermeidate ouputs used by germline pipeline (such as the output BAM) otherwise they won't be kept as outputs."
  }

  String pipeline_version = "3.2.3"

  References references = alignment_references.references


  call InternalTasks.MakeSafeFilename as MakeSafeFilename {
    input:
      name = base_file_name
  }

  call Tasks.VerifyPipelineInputs as VerifyPipelineInputs {
    input:
      input_cram_list = input_cram_list,
      input_bam_list  = input_bam_list
  }

  call UltimaGenomicsWholeGenomeGermlineAlignmentMarkDuplicates.AlignmentAndMarkDuplicates as AlignmentAndMarkDuplicates {
    input:
      input_cram_bam                = select_first([input_cram_list,input_bam_list]),
      is_cram                       = VerifyPipelineInputs.is_cram,
      base_file_name_sub            = MakeSafeFilename.output_safe_name,
      reads_per_split               = reads_per_split,
      rsq_threshold                 = rsq_threshold,
      alignment_references          = alignment_references,
      references                    = references
  }

  # Convert the final merged recalibrated BAM file to CRAM format
  call Utilities.ConvertToCram {
    input:
      input_bam       = AlignmentAndMarkDuplicates.output_bam,
      ref_fasta       = references.ref_fasta,
      ref_fasta_index = references.ref_fasta_index,
      output_basename = MakeSafeFilename.output_safe_name
  }

  Float dynamic_validate_cram_disk_size = size(ConvertToCram.output_cram, "GB") + size(references.ref_fasta, "GB") + size(references.ref_fasta_index, "GB") + size(references.ref_dict, "GB") + 80
  Int validate_cram_disk_size = if dynamic_validate_cram_disk_size > 510 then ceil(dynamic_validate_cram_disk_size) else 510

  # Validate the CRAM file
  call QC.ValidateSamFile as ValidateCram {
    input:
      input_bam         = ConvertToCram.output_cram,
      input_bam_index   = ConvertToCram.output_cram_index,
      report_filename   = MakeSafeFilename.output_safe_name + ".cram.validation_report",
      ref_dict          = references.ref_dict,
      ref_fasta         = references.ref_fasta,
      ref_fasta_index   = references.ref_fasta_index,
      ignore            = ["MISSING_TAG_NM" ,"INVALID_PLATFORM_VALUE"],
      max_output        = 1000000000,
      is_outlier_data   = true, #sets SKIP_MATE_VALIDATION=true
      disk_size         = validate_cram_disk_size
  }

  call Tasks.ExtractSampleNameFlowOrder {
    input:
      input_bam  = AlignmentAndMarkDuplicates.output_bam,
      references = references
  }

  call UltimaGenomicsWholeGenomeGermlineQC.UltimaGenomicsWholeGenomeGermlineQC as CollectStatistics {
    input:
      agg_bam                               = AlignmentAndMarkDuplicates.output_bam,
      agg_bam_index                         = AlignmentAndMarkDuplicates.output_bam_index,
      base_file_name                        = MakeSafeFilename.output_safe_name,
      base_file_name_sub                    = MakeSafeFilename.output_safe_name,
      references                            = references,
      contamination_sites                   = contamination_sites,
      wgs_coverage_interval_list            = vcf_post_processing.wgs_coverage_interval_list,
      max_duplication_in_reasonable_sample  = vcf_post_processing.max_duplication_in_reasonable_sample,
      max_chimerism_in_reasonable_sample    = vcf_post_processing.max_chimerism_in_reasonable_sample,
      flow_order                            = ExtractSampleNameFlowOrder.flow_order
  }

  call Utilities.MakeOptionalOutputBam {
    input:
      bam_input = AlignmentAndMarkDuplicates.output_bam,
      bai_input = AlignmentAndMarkDuplicates.output_bam_index,
      keep_inputs = save_bam_file
  }

  # Outputs that will be retained when execution is complete
  output {
    File output_cram = ConvertToCram.output_cram
    File output_cram_index = ConvertToCram.output_cram_index
    File output_cram_md5 = ConvertToCram.output_cram_md5

    File selfSM = CollectStatistics.selfSM
    Float contamination = CollectStatistics.contamination

    # STATISTIC COLLECTION
    File quality_yield_metrics = CollectStatistics.quality_yield_metrics
    File wgs_metrics = CollectStatistics.wgs_metrics
    File raw_wgs_metrics = CollectStatistics.raw_wgs_metrics
    File duplicate_metrics = CollectStatistics.duplicate_metrics
    File agg_alignment_summary_metrics = CollectStatistics.agg_alignment_summary_metrics
    File? agg_alignment_summary_pdf = CollectStatistics.agg_alignment_summary_pdf
    File agg_gc_bias_detail_metrics = CollectStatistics.agg_gc_bias_detail_metrics
    File agg_gc_bias_pdf = CollectStatistics.agg_gc_bias_pdf
    File agg_gc_bias_summary_metrics = CollectStatistics.agg_gc_bias_summary_metrics
    File agg_quality_distribution_pdf = CollectStatistics.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = CollectStatistics.agg_quality_distribution_metrics
    Float duplication_rate = CollectStatistics.duplication_rate
    Float chimerism_rate = CollectStatistics.chimerism_rate
    Boolean is_outlier_data = CollectStatistics.is_outlier_data

    String sample_name = ExtractSampleNameFlowOrder.sample_name
    String flow_order = ExtractSampleNameFlowOrder.flow_order
    String barcode = ExtractSampleNameFlowOrder.barcode_seq
    String id = ExtractSampleNameFlowOrder.readgroup_id

    #Intermediate outputs required for germline pipeline
    File? output_bam = MakeOptionalOutputBam.optional_output_bam
    File? output_bam_index = MakeOptionalOutputBam.optional_output_bai

    String output_safe_name = MakeSafeFilename.output_safe_name
  }

}
