version 1.0

import "../../../../../../tasks/wdl/UltimaGenomicsWholeGenomeGermlineTasks.wdl" as Tasks
import "../../../../../../tasks/wdl/Utilities.wdl" as Utilities
import "../../../../../../tasks/wdl/GermlineVariantDiscovery.wdl" as VariantDiscoverTasks
import "../../../../../../tasks/wdl/UltimaGenomicsWholeGenomeGermlineAlignmentMarkDuplicates.wdl" as UltimaGenomicsWholeGenomeGermlineAlignmentMarkDuplicates
import "../../../../../../tasks/wdl/InternalTasks.wdl" as InternalTasks
import "../../../../../../tasks/wdl/Qc.wdl" as QC
import "../../../../../../tasks/wdl/UltimaGenomicsWholeGenomeGermlineQC.wdl" as UltimaGenomicsWholeGenomeGermlineQC
import "../../../../../../structs/dna_seq/UltimaGenomicsWholeGenomeGermlineStructs.wdl" as Structs


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

  String pipeline_version = "1.1.3"

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
      references                    = references,
      save_bam_file                 = save_bam_file
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
    File? output_bam = AlignmentAndMarkDuplicates.optional_output_bam
    File? output_bam_index = AlignmentAndMarkDuplicates.optional_output_bam_index

    String output_safe_name = MakeSafeFilename.output_safe_name
  }

}
