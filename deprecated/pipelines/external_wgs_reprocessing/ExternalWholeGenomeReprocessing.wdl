version 1.0
# ExternalWholeGenomeReprocessing is now deprecated
import "../../../../../pipelines/wdl/reprocessing/wgs/WholeGenomeReprocessing.wdl" as WholeGenomeReprocessing
import "../../../../../tasks/wdl/CopyFilesFromCloudToCloud.wdl" as Copy

workflow ExternalWholeGenomeReprocessing {


  String pipeline_version = "2.3.4"

  input {
    File? input_cram
    File? input_bam
    File? output_map

    String sample_name
    String base_file_name
    String final_gvcf_base_name
    String unmapped_bam_suffix

    File cram_ref_fasta
    File cram_ref_fasta_index

    File wgs_coverage_interval_list

    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index

    DNASeqSingleSampleReferences references
    VariantCallingScatterSettings scatter_settings
    PapiSettings papi_settings

    String destination_cloud_path
    String vault_token_path
    String google_account_vault_path

    String cloud_provider
  }

  call WholeGenomeReprocessing.WholeGenomeReprocessing {
    input:
    input_cram = input_cram,
    input_bam = input_bam,
    output_map = output_map,
    sample_name = sample_name,
    base_file_name = base_file_name,
    final_gvcf_base_name = final_gvcf_base_name,
    unmapped_bam_suffix = unmapped_bam_suffix,
    cram_ref_fasta = cram_ref_fasta,
    cram_ref_fasta_index = cram_ref_fasta_index,
    references = references,
    fingerprint_genotypes_file = fingerprint_genotypes_file,
    fingerprint_genotypes_index = fingerprint_genotypes_index,
    papi_settings = papi_settings,
    wgs_coverage_interval_list = wgs_coverage_interval_list,
    scatter_settings = scatter_settings,
    cloud_provider = cloud_provider
  }

  call Copy.CopyFilesFromCloudToCloud {
    input:
      files_to_copy = flatten([
                        # The File outputs
                        [WholeGenomeReprocessing.read_group_alignment_summary_metrics,
                          WholeGenomeReprocessing.selfSM,
                          WholeGenomeReprocessing.calculate_read_group_checksum_md5,
                          WholeGenomeReprocessing.agg_alignment_summary_metrics,
                          WholeGenomeReprocessing.agg_bait_bias_detail_metrics,
                          WholeGenomeReprocessing.agg_bait_bias_summary_metrics,
                          WholeGenomeReprocessing.agg_gc_bias_detail_metrics,
                          WholeGenomeReprocessing.agg_gc_bias_pdf,
                          WholeGenomeReprocessing.agg_gc_bias_summary_metrics,
                          WholeGenomeReprocessing.agg_insert_size_histogram_pdf,
                          WholeGenomeReprocessing.agg_insert_size_metrics,
                          WholeGenomeReprocessing.agg_pre_adapter_detail_metrics,
                          WholeGenomeReprocessing.agg_pre_adapter_summary_metrics,
                          WholeGenomeReprocessing.agg_quality_distribution_pdf,
                          WholeGenomeReprocessing.agg_quality_distribution_metrics,
                          WholeGenomeReprocessing.duplicate_metrics,
                          WholeGenomeReprocessing.gvcf_summary_metrics,
                          WholeGenomeReprocessing.gvcf_detail_metrics,
                          WholeGenomeReprocessing.wgs_metrics,
                          WholeGenomeReprocessing.raw_wgs_metrics,
                          WholeGenomeReprocessing.output_cram,
                          WholeGenomeReprocessing.output_cram_index,
                          WholeGenomeReprocessing.output_cram_md5,
                          WholeGenomeReprocessing.validate_cram_file_report,
                          WholeGenomeReprocessing.output_vcf,
                          WholeGenomeReprocessing.output_vcf_index],
                        # The Array[File] outputs
                        WholeGenomeReprocessing.validation_report,
                        WholeGenomeReprocessing.quality_yield_metrics,
                        WholeGenomeReprocessing.unsorted_read_group_base_distribution_by_cycle_pdf,
                        WholeGenomeReprocessing.unsorted_read_group_base_distribution_by_cycle_metrics,
                        WholeGenomeReprocessing.unsorted_read_group_insert_size_histogram_pdf,
                        WholeGenomeReprocessing.unsorted_read_group_insert_size_metrics,
                        WholeGenomeReprocessing.unsorted_read_group_quality_by_cycle_pdf,
                        WholeGenomeReprocessing.unsorted_read_group_quality_by_cycle_metrics,
                        WholeGenomeReprocessing.unsorted_read_group_quality_distribution_pdf,
                        WholeGenomeReprocessing.unsorted_read_group_quality_distribution_metrics,
                        # The File? outputs
                        select_all([WholeGenomeReprocessing.cross_check_fingerprints_metrics]),
                        select_all([WholeGenomeReprocessing.fingerprint_summary_metrics]),
                        select_all([WholeGenomeReprocessing.fingerprint_detail_metrics]),
                        select_all([WholeGenomeReprocessing.output_bqsr_reports])]),
        vault_token_path = vault_token_path,
        destination_cloud_path = destination_cloud_path,
        google_account_vault_path = google_account_vault_path,
        contamination = WholeGenomeReprocessing.contamination,
        base_file_name = base_file_name
  }

  output {
    Array[File] validation_report = WholeGenomeReprocessing.validation_report
    Array[File] unmapped_bams = WholeGenomeReprocessing.unmapped_bams

    Array[File] quality_yield_metrics = WholeGenomeReprocessing.quality_yield_metrics

    Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = WholeGenomeReprocessing.unsorted_read_group_base_distribution_by_cycle_pdf
    Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = WholeGenomeReprocessing.unsorted_read_group_base_distribution_by_cycle_metrics
    Array[File] unsorted_read_group_insert_size_histogram_pdf = WholeGenomeReprocessing.unsorted_read_group_insert_size_histogram_pdf
    Array[File] unsorted_read_group_insert_size_metrics = WholeGenomeReprocessing.unsorted_read_group_insert_size_metrics
    Array[File] unsorted_read_group_quality_by_cycle_pdf = WholeGenomeReprocessing.unsorted_read_group_quality_by_cycle_pdf
    Array[File] unsorted_read_group_quality_by_cycle_metrics = WholeGenomeReprocessing.unsorted_read_group_quality_by_cycle_metrics
    Array[File] unsorted_read_group_quality_distribution_pdf = WholeGenomeReprocessing.unsorted_read_group_quality_distribution_pdf
    Array[File] unsorted_read_group_quality_distribution_metrics = WholeGenomeReprocessing.unsorted_read_group_quality_distribution_metrics

    File read_group_alignment_summary_metrics = WholeGenomeReprocessing.read_group_alignment_summary_metrics
    File read_group_gc_bias_detail_metrics = WholeGenomeReprocessing.read_group_gc_bias_detail_metrics
    File read_group_gc_bias_pdf = WholeGenomeReprocessing.read_group_gc_bias_pdf
    File read_group_gc_bias_summary_metrics = WholeGenomeReprocessing.read_group_gc_bias_summary_metrics

    File? cross_check_fingerprints_metrics = WholeGenomeReprocessing.cross_check_fingerprints_metrics

    File selfSM = WholeGenomeReprocessing.selfSM
    Float contamination = WholeGenomeReprocessing.contamination

    File calculate_read_group_checksum_md5 = WholeGenomeReprocessing.calculate_read_group_checksum_md5

    File agg_alignment_summary_metrics = WholeGenomeReprocessing.agg_alignment_summary_metrics
    File agg_bait_bias_detail_metrics = WholeGenomeReprocessing.agg_bait_bias_detail_metrics
    File agg_bait_bias_summary_metrics = WholeGenomeReprocessing.agg_bait_bias_summary_metrics
    File agg_gc_bias_detail_metrics = WholeGenomeReprocessing.agg_gc_bias_detail_metrics
    File agg_gc_bias_pdf = WholeGenomeReprocessing.agg_gc_bias_pdf
    File agg_gc_bias_summary_metrics = WholeGenomeReprocessing.agg_gc_bias_summary_metrics
    File agg_insert_size_histogram_pdf = WholeGenomeReprocessing.agg_insert_size_histogram_pdf
    File agg_insert_size_metrics = WholeGenomeReprocessing.agg_insert_size_metrics
    File agg_pre_adapter_detail_metrics = WholeGenomeReprocessing.agg_pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = WholeGenomeReprocessing.agg_pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = WholeGenomeReprocessing.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = WholeGenomeReprocessing.agg_quality_distribution_metrics

    File? fingerprint_summary_metrics = WholeGenomeReprocessing.fingerprint_summary_metrics
    File? fingerprint_detail_metrics = WholeGenomeReprocessing.fingerprint_detail_metrics

    File wgs_metrics = WholeGenomeReprocessing.wgs_metrics
    File raw_wgs_metrics = WholeGenomeReprocessing.raw_wgs_metrics

    File duplicate_metrics = WholeGenomeReprocessing.duplicate_metrics
    File? output_bqsr_reports = WholeGenomeReprocessing.output_bqsr_reports

    File gvcf_summary_metrics = WholeGenomeReprocessing.gvcf_summary_metrics
    File gvcf_detail_metrics = WholeGenomeReprocessing.gvcf_detail_metrics

    File output_cram = WholeGenomeReprocessing.output_cram
    File output_cram_index = WholeGenomeReprocessing.output_cram_index
    File output_cram_md5 = WholeGenomeReprocessing.output_cram_md5

    File validate_cram_file_report = WholeGenomeReprocessing.validate_cram_file_report

    File output_vcf = WholeGenomeReprocessing.output_vcf
    File output_vcf_index = WholeGenomeReprocessing.output_vcf_index
  }
  meta {
    allowNestedInputs: true
  }
}
