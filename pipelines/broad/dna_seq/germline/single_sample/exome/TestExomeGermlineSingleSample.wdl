version 1.0

import "../../../../../../pipelines/broad/dna_seq/germline/single_sample/exome/ExomeGermlineSingleSample.wdl" as ExomeGermlineSingleSample
import "../../../../../../verification/VerifyGermlineSingleSample.wdl" as VerifyGermlineSingleSample
import "../../../../../../tasks/broad/Utilities.wdl" as Utilities
import "../../../../../../tasks/broad/CopyFilesFromCloudToCloud.wdl" as Copy

workflow TestExomeGermlineSingleSample {

  input {
    PapiSettings papi_settings
    SampleAndUnmappedBams sample_and_unmapped_bams
    DNASeqSingleSampleReferences references
    VariantCallingScatterSettings scatter_settings

    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index

    File target_interval_list
    File bait_interval_list
    String bait_set_name

    Boolean provide_bam_output = false

    # These values will be determined and injected into the inputs by the scala test framework
    String? truth_path
    String results_path
    Boolean? use_timestamp
    Boolean? update_truth
    String? timestamp
    String cromwell_url_auth
    String vault_token_path
    String google_account_vault_path
    #Array[String] metrics_files_to_test
    #Boolean update_truth
  }

  meta {
    allowNestedInputs: true
  }

  # Run the pipeline
  call ExomeGermlineSingleSample.ExomeGermlineSingleSample {
    input:
      sample_and_unmapped_bams     = sample_and_unmapped_bams,
      references                   = references,
      scatter_settings             = scatter_settings,
      fingerprint_genotypes_file   = fingerprint_genotypes_file,
      fingerprint_genotypes_index  = fingerprint_genotypes_index,
      papi_settings                = papi_settings,
      target_interval_list         = target_interval_list,
      bait_interval_list           = bait_interval_list,
      bait_set_name                = bait_set_name,
  }

  call Copy.CopyFilesFromCloudToCloud as CopyToTestResults {
    input:
      files_to_copy = flatten([
                            [ # File outputs
                            ExomeGermlineSingleSample.read_group_alignment_summary_metrics,
                            ExomeGermlineSingleSample.selfSM,
                            ExomeGermlineSingleSample.calculate_read_group_checksum_md5,
                            ExomeGermlineSingleSample.agg_alignment_summary_metrics,
                            ExomeGermlineSingleSample.agg_bait_bias_detail_metrics,
                            ExomeGermlineSingleSample.agg_bait_bias_summary_metrics,
                            ExomeGermlineSingleSample.agg_insert_size_histogram_pdf,
                            ExomeGermlineSingleSample.agg_insert_size_metrics,
                            ExomeGermlineSingleSample.agg_pre_adapter_detail_metrics,
                            ExomeGermlineSingleSample.agg_pre_adapter_summary_metrics,
                            ExomeGermlineSingleSample.agg_quality_distribution_pdf,
                            ExomeGermlineSingleSample.agg_quality_distribution_metrics,
                            ExomeGermlineSingleSample.agg_error_summary_metrics,
                            ExomeGermlineSingleSample.duplicate_metrics,
                            ExomeGermlineSingleSample.gvcf_summary_metrics,
                            ExomeGermlineSingleSample.gvcf_detail_metrics,
                            ExomeGermlineSingleSample.hybrid_selection_metrics,
                            ExomeGermlineSingleSample.output_cram,
                            ExomeGermlineSingleSample.output_cram_index,
                            ExomeGermlineSingleSample.output_cram_md5,
                            ExomeGermlineSingleSample.validate_cram_file_report,
                            ExomeGermlineSingleSample.output_vcf,
                            ExomeGermlineSingleSample.output_vcf_index
                            ], # Array[File] outputs
                            ExomeGermlineSingleSample.quality_yield_metrics,
                            ExomeGermlineSingleSample.unsorted_read_group_base_distribution_by_cycle_pdf,
                            ExomeGermlineSingleSample.unsorted_read_group_base_distribution_by_cycle_metrics,
                            ExomeGermlineSingleSample.unsorted_read_group_insert_size_histogram_pdf,
                            ExomeGermlineSingleSample.unsorted_read_group_insert_size_metrics,
                            ExomeGermlineSingleSample.unsorted_read_group_quality_by_cycle_pdf,
                            ExomeGermlineSingleSample.unsorted_read_group_quality_by_cycle_metrics,
                            ExomeGermlineSingleSample.unsorted_read_group_quality_distribution_pdf,
                            ExomeGermlineSingleSample.unsorted_read_group_quality_distribution_metrics,
                            # File? outputs
                            select_all([ExomeGermlineSingleSample.cross_check_fingerprints_metrics]),
                            select_all([ExomeGermlineSingleSample.fingerprint_summary_metrics]),
                            select_all([ExomeGermlineSingleSample.fingerprint_detail_metrics]),
                            select_all([ExomeGermlineSingleSample.output_bqsr_reports]),
                            select_all([ExomeGermlineSingleSample.output_bam]),
                            select_all([ExomeGermlineSingleSample.output_bam_index]),
      ]),
      vault_token_path = vault_token_path,
      google_account_vault_path = google_account_vault_path,
      contamination = ExomeGermlineSingleSample.contamination,
      destination_cloud_path = results_path
  }

  ## Copy the outputs from broad-gotc-cromwell-execution
  ## to broad-gotc-test-results
  #call Utilities.CopyWorkflowOutputsByPath as CopyToTestResults {
  #  input:
  #    output_file_path          = ExomeGermlineSingleSample.output_vcf,
  #    copy_bucket_path          = results_path,
  #    workflow_name             = "ExomeGermlineSingleSample",
  #    cromwell_url              = cromwell_url_auth,
  #    vault_token_path          = vault_token_path,
  #    google_account_vault_path = google_account_vault_path
#
  #}


}
  # Outputs are intentionally left undefined to automatically propagate all outputs from all tasks (without needing to individually define the outputs)
