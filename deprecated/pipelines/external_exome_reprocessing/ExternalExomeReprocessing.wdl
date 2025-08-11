version 1.0

# ExternalExomeReprocessing is now deprecated 2025-03-06

import "../../../../../pipelines/wdl/reprocessing/exome/ExomeReprocessing.wdl" as ExomeReprocessing
import "../../../../../tasks/wdl/CopyFilesFromCloudToCloud.wdl" as Copy

workflow ExternalExomeReprocessing {

  String pipeline_version = "3.3.4"


  input {
    File? input_cram
    File? input_bam

    String sample_name
    String base_file_name
    String final_gvcf_base_name
    String unmapped_bam_suffix

    File cram_ref_fasta
    File cram_ref_fasta_index

    String bait_set_name
    File bait_interval_list
    File target_interval_list

    DNASeqSingleSampleReferences references
    VariantCallingScatterSettings scatter_settings
    PapiSettings papi_settings

    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index

    String destination_cloud_path
    String vault_token_path
    String google_account_vault_path

    String cloud_provider
  }

  call ExomeReprocessing.ExomeReprocessing {
    input:
      input_cram = input_cram,
      input_bam = input_bam,
      sample_name = sample_name,
      base_file_name = base_file_name,
      final_gvcf_base_name = final_gvcf_base_name,
      unmapped_bam_suffix = unmapped_bam_suffix,
      bait_set_name = bait_set_name,
      bait_interval_list = bait_interval_list,
      target_interval_list = target_interval_list,
      references = references,
      scatter_settings = scatter_settings,
      fingerprint_genotypes_file = fingerprint_genotypes_file,
      fingerprint_genotypes_index = fingerprint_genotypes_index,
      cram_ref_fasta = cram_ref_fasta,
      cram_ref_fasta_index = cram_ref_fasta_index,
      papi_settings = papi_settings,
      cloud_provider = cloud_provider
  }

  call Copy.CopyFilesFromCloudToCloud {
    input:
      files_to_copy = flatten([
                        # The File outputs
                        [ExomeReprocessing.read_group_alignment_summary_metrics,
                          ExomeReprocessing.selfSM,
                          ExomeReprocessing.calculate_read_group_checksum_md5,
                          ExomeReprocessing.agg_alignment_summary_metrics,
                          ExomeReprocessing.agg_bait_bias_detail_metrics,
                          ExomeReprocessing.agg_bait_bias_summary_metrics,
                          ExomeReprocessing.agg_insert_size_histogram_pdf,
                          ExomeReprocessing.agg_insert_size_metrics,
                          ExomeReprocessing.agg_pre_adapter_detail_metrics,
                          ExomeReprocessing.agg_pre_adapter_summary_metrics,
                          ExomeReprocessing.agg_quality_distribution_pdf,
                          ExomeReprocessing.agg_quality_distribution_metrics,
                          ExomeReprocessing.duplicate_metrics,
                          ExomeReprocessing.gvcf_summary_metrics,
                          ExomeReprocessing.gvcf_detail_metrics,
                          ExomeReprocessing.hybrid_selection_metrics,
                          ExomeReprocessing.output_cram,
                          ExomeReprocessing.output_cram_index,
                          ExomeReprocessing.output_cram_md5,
                          ExomeReprocessing.validate_cram_file_report,
                          ExomeReprocessing.output_vcf,
                          ExomeReprocessing.output_vcf_index],
                        # The Array[File] outputs
                        ExomeReprocessing.validation_report,
                        ExomeReprocessing.quality_yield_metrics,
                        ExomeReprocessing.unsorted_read_group_base_distribution_by_cycle_pdf,
                        ExomeReprocessing.unsorted_read_group_base_distribution_by_cycle_metrics,
                        ExomeReprocessing.unsorted_read_group_insert_size_histogram_pdf,
                        ExomeReprocessing.unsorted_read_group_insert_size_metrics,
                        ExomeReprocessing.unsorted_read_group_quality_by_cycle_pdf,
                        ExomeReprocessing.unsorted_read_group_quality_by_cycle_metrics,
                        ExomeReprocessing.unsorted_read_group_quality_distribution_pdf,
                        ExomeReprocessing.unsorted_read_group_quality_distribution_metrics,
                        # The File? outputs
                        select_all([ExomeReprocessing.cross_check_fingerprints_metrics]),
                        select_all([ExomeReprocessing.fingerprint_summary_metrics]),
                        select_all([ExomeReprocessing.fingerprint_detail_metrics]),
                        select_all([ExomeReprocessing.output_bqsr_reports])]),
        vault_token_path = vault_token_path,
        destination_cloud_path = destination_cloud_path,
        google_account_vault_path = google_account_vault_path,
        contamination = ExomeReprocessing.contamination,
        base_file_name = base_file_name
  }

  output {
    Array[File] sam_validation_report = ExomeReprocessing.validation_report

    Array[File] quality_yield_metrics = ExomeReprocessing.quality_yield_metrics

    Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = ExomeReprocessing.unsorted_read_group_base_distribution_by_cycle_pdf
    Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = ExomeReprocessing.unsorted_read_group_base_distribution_by_cycle_metrics
    Array[File] unsorted_read_group_insert_size_histogram_pdf = ExomeReprocessing.unsorted_read_group_insert_size_histogram_pdf
    Array[File] unsorted_read_group_insert_size_metrics = ExomeReprocessing.unsorted_read_group_insert_size_metrics
    Array[File] unsorted_read_group_quality_by_cycle_pdf = ExomeReprocessing.unsorted_read_group_quality_by_cycle_pdf
    Array[File] unsorted_read_group_quality_by_cycle_metrics = ExomeReprocessing.unsorted_read_group_quality_by_cycle_metrics
    Array[File] unsorted_read_group_quality_distribution_pdf = ExomeReprocessing.unsorted_read_group_quality_distribution_pdf
    Array[File] unsorted_read_group_quality_distribution_metrics = ExomeReprocessing.unsorted_read_group_quality_distribution_metrics

    File read_group_alignment_summary_metrics = ExomeReprocessing.read_group_alignment_summary_metrics

    File? cross_check_fingerprints_metrics = ExomeReprocessing.cross_check_fingerprints_metrics

    File selfSM = ExomeReprocessing.selfSM
    Float contamination = ExomeReprocessing.contamination

    File calculate_read_group_checksum_md5 = ExomeReprocessing.calculate_read_group_checksum_md5

    File agg_alignment_summary_metrics = ExomeReprocessing.agg_alignment_summary_metrics
    File agg_bait_bias_detail_metrics = ExomeReprocessing.agg_bait_bias_detail_metrics
    File agg_bait_bias_summary_metrics = ExomeReprocessing.agg_bait_bias_summary_metrics
    File agg_insert_size_histogram_pdf = ExomeReprocessing.agg_insert_size_histogram_pdf
    File agg_insert_size_metrics = ExomeReprocessing.agg_insert_size_metrics
    File agg_pre_adapter_detail_metrics = ExomeReprocessing.agg_pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = ExomeReprocessing.agg_pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = ExomeReprocessing.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = ExomeReprocessing.agg_quality_distribution_metrics

    File? fingerprint_summary_metrics = ExomeReprocessing.fingerprint_summary_metrics
    File? fingerprint_detail_metrics = ExomeReprocessing.fingerprint_detail_metrics

    File duplicate_metrics = ExomeReprocessing.duplicate_metrics
    File? output_bqsr_reports = ExomeReprocessing.output_bqsr_reports

    File gvcf_summary_metrics = ExomeReprocessing.gvcf_summary_metrics
    File gvcf_detail_metrics = ExomeReprocessing.gvcf_detail_metrics

    File hybrid_selection_metrics = ExomeReprocessing.hybrid_selection_metrics

    File output_cram = ExomeReprocessing.output_cram
    File output_cram_index = ExomeReprocessing.output_cram_index
    File output_cram_md5 = ExomeReprocessing.output_cram_md5

    File validate_cram_file_report = ExomeReprocessing.validate_cram_file_report

    File output_vcf = ExomeReprocessing.output_vcf
    File output_vcf_index = ExomeReprocessing.output_vcf_index
  }
  meta {
    allowNestedInputs: true
  }
}
