version 1.0

import "../../pipelines/broad/dna_seq/germline/single_sample/exome/ExomeGermlineSingleSample.wdl" as ExomeGermlineSingleSample
import "../../verification/VerifyGermlineSingleSample.wdl" as VerifyGermlineSingleSample
import "../../tasks/broad/Utilities.wdl" as Utilities
import "../../tasks/broad/CopyFilesFromCloudToCloud.wdl" as Copy

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
    String truth_path
    String results_path
    Boolean update_truth
    String vault_token_path
    String google_account_vault_path
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
      provide_bam_output           = provide_bam_output
  }

  # Collect all of the pipeline outputs into a single Array[String]]
  Array[String] pipeline_outputs = flatten([
                            [ # File outputs
                            ExomeGermlineSingleSample.selfSM,
                            ExomeGermlineSingleSample.agg_insert_size_histogram_pdf,
                            ExomeGermlineSingleSample.agg_quality_distribution_pdf,
                            ExomeGermlineSingleSample.calculate_read_group_checksum_md5,
                            ExomeGermlineSingleSample.agg_insert_size_histogram_pdf,
                            ExomeGermlineSingleSample.agg_quality_distribution_pdf,
                            ExomeGermlineSingleSample.output_cram,
                            ExomeGermlineSingleSample.output_cram_index,
                            ExomeGermlineSingleSample.output_cram_md5,
                            ExomeGermlineSingleSample.validate_cram_file_report
                            ], # Array[File] outputs
                            ExomeGermlineSingleSample.unsorted_read_group_base_distribution_by_cycle_pdf,
                            ExomeGermlineSingleSample.unsorted_read_group_insert_size_histogram_pdf,
                            ExomeGermlineSingleSample.unsorted_read_group_quality_by_cycle_pdf,
                            ExomeGermlineSingleSample.unsorted_read_group_quality_distribution_pdf,
                            # File? outputs
                            select_all([ExomeGermlineSingleSample.output_bqsr_reports]),
                            select_all([ExomeGermlineSingleSample.output_bam]),
                            select_all([ExomeGermlineSingleSample.output_bam_index]),
  ])

  # Collect all of the pipeline metrics into a single Array[String]
  Array[String] pipeline_metrics = flatten([
                              [ # File outputs
                              ExomeGermlineSingleSample.read_group_alignment_summary_metrics,
                              ExomeGermlineSingleSample.agg_alignment_summary_metrics,
                              ExomeGermlineSingleSample.agg_bait_bias_detail_metrics,
                              ExomeGermlineSingleSample.agg_bait_bias_summary_metrics,
                              ExomeGermlineSingleSample.agg_insert_size_metrics,
                              ExomeGermlineSingleSample.agg_pre_adapter_detail_metrics,
                              ExomeGermlineSingleSample.agg_pre_adapter_summary_metrics,
                              ExomeGermlineSingleSample.agg_quality_distribution_metrics,
                              ExomeGermlineSingleSample.agg_error_summary_metrics,
                              ExomeGermlineSingleSample.duplicate_metrics,
                              ExomeGermlineSingleSample.hybrid_selection_metrics,
                              ], # Array[File] outputs
                              ExomeGermlineSingleSample.quality_yield_metrics,
                              ExomeGermlineSingleSample.unsorted_read_group_base_distribution_by_cycle_metrics,
                              ExomeGermlineSingleSample.unsorted_read_group_insert_size_metrics,
                              ExomeGermlineSingleSample.unsorted_read_group_quality_by_cycle_metrics,
                              ExomeGermlineSingleSample.unsorted_read_group_quality_distribution_metrics,
                              # File? outputs
                              select_all([ExomeGermlineSingleSample.cross_check_fingerprints_metrics]),
                              select_all([ExomeGermlineSingleSample.fingerprint_summary_metrics]),
                              select_all([ExomeGermlineSingleSample.fingerprint_detail_metrics]),
  ])

  # Copy results of pipeline to test results bucket
  call Copy.CopyFilesFromCloudToCloud as CopyToTestResults {
    input:
      files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
      vault_token_path          = vault_token_path,
      google_account_vault_path = google_account_vault_path,
      contamination             = ExomeGermlineSingleSample.contamination,
      destination_cloud_path    = results_path
  }

  # If updating truth then copy pipeline results to truth bucket
  if (update_truth){
    call Copy.CopyFilesFromCloudToCloud as CopyToTruth {
      input:
        files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
        vault_token_path          = vault_token_path,
        google_account_vault_path = google_account_vault_path,
        contamination             = ExomeGermlineSingleSample.contamination,
        destination_cloud_path    = truth_path
      }
  }

  # If not updating truth then we need to collect all input for the validation WDL
  # This is achieved by passing each desired file/array[files] to GetValidationInputs
  if (!update_truth){
    call Utilities.GetValidationInputs as GetMetricsInputs {
      input:
        input_files  = pipeline_metrics,
        results_path = results_path,
        truth_path   = truth_path
    }

    call Utilities.GetValidationInputs as GetCrams {
      input:
          input_file   = ExomeGermlineSingleSample.output_cram,
          results_path = results_path,
          truth_path   = truth_path
    }

    call Utilities.GetValidationInputs as GetCrais {
      input:
        input_file   = ExomeGermlineSingleSample.output_cram_index,
        results_path = results_path,
        truth_path   = truth_path
    }
    
    # done is dummy input to force copy completion before verification
    call VerifyGermlineSingleSample.VerifyGermlineSingleSample as Verify {
      input:
        truth_metrics    = GetMetricsInputs.truth_files,
        truth_cram       = GetCrams.truth_file,
        truth_crai       = GetCrais.truth_file,
        test_metrics     = GetMetricsInputs.results_files,
        test_cram        = GetCrams.results_file,
        test_crai        = GetCrais.results_file,
        done             = CopyToTestResults.done
    }
  }

  output {
    Array[File]? metric_comparison_report_files = Verify.metric_comparison_report_files
  }

}
