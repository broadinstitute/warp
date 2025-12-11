version 1.0


import "../../pipelines/wdl/reprocessing/external/exome/ExternalExomeReprocessing.wdl" as ExternalExomeReprocessing
import "../../verification/VerifyExternalReprocessing.wdl" as VerifyExternalReprocessing
import "../../tasks/wdl/CopyFilesFromCloudToCloud.wdl" as Copy

workflow TestExternalExomeReprocessing {

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
      # External pipelines copy to a desired cloud path after running reprocessing
      String destination_cloud_path
      String cloud_provider

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
  
    call ExternalExomeReprocessing.ExternalExomeReprocessing {
      input:
        input_cram = input_cram,
        input_bam = input_bam,
        sample_name = sample_name,
        base_file_name = base_file_name,
        final_gvcf_base_name = final_gvcf_base_name,
        unmapped_bam_suffix = unmapped_bam_suffix,
        cram_ref_fasta = cram_ref_fasta,
        cram_ref_fasta_index = cram_ref_fasta_index,
        bait_set_name = bait_set_name,
        bait_interval_list = bait_interval_list,
        target_interval_list = target_interval_list,
        references = references,
        scatter_settings = scatter_settings,
        papi_settings = papi_settings,
        fingerprint_genotypes_file = fingerprint_genotypes_file,
        fingerprint_genotypes_index = fingerprint_genotypes_index,
        destination_cloud_path = destination_cloud_path,
        vault_token_path = vault_token_path,
        google_account_vault_path = google_account_vault_path,
        cloud_provider = cloud_provider
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    ExternalExomeReprocessing.output_vcf_index,
                                    ExternalExomeReprocessing.output_vcf,
                                    ExternalExomeReprocessing.validate_cram_file_report,
                                    ExternalExomeReprocessing.output_cram_md5,
                                    ExternalExomeReprocessing.output_cram_index,
                                    ExternalExomeReprocessing.output_cram,
                                    ExternalExomeReprocessing.agg_quality_distribution_pdf,
                                    ExternalExomeReprocessing.agg_insert_size_histogram_pdf,
                                    ExternalExomeReprocessing.calculate_read_group_checksum_md5,
                                    ExternalExomeReprocessing.selfSM,
                                    ],
                                    # Array[File] outputs
                                    ExternalExomeReprocessing.unsorted_read_group_quality_distribution_pdf,
                                    ExternalExomeReprocessing.unsorted_read_group_quality_by_cycle_pdf,
                                    ExternalExomeReprocessing.unsorted_read_group_insert_size_histogram_pdf,
                                    ExternalExomeReprocessing.unsorted_read_group_base_distribution_by_cycle_pdf,
                                    ExternalExomeReprocessing.sam_validation_report,
                                    # File? outputs
                                    select_all([ExternalExomeReprocessing.output_bqsr_reports]),
                                    
    ])

    
    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    [ # File outputs
                                    ExternalExomeReprocessing.hybrid_selection_metrics,
                                    ExternalExomeReprocessing.gvcf_detail_metrics,
                                    ExternalExomeReprocessing.gvcf_summary_metrics,
                                    ExternalExomeReprocessing.duplicate_metrics,
                                    ExternalExomeReprocessing.agg_quality_distribution_metrics,
                                    ExternalExomeReprocessing.agg_pre_adapter_summary_metrics,
                                    ExternalExomeReprocessing.agg_pre_adapter_detail_metrics,
                                    ExternalExomeReprocessing.agg_insert_size_metrics,
                                    ExternalExomeReprocessing.agg_bait_bias_summary_metrics,
                                    ExternalExomeReprocessing.agg_bait_bias_detail_metrics,
                                    ExternalExomeReprocessing.agg_alignment_summary_metrics,
                                    ExternalExomeReprocessing.read_group_alignment_summary_metrics,
                                    ],
                                    # Array[File] outputs
                                    ExternalExomeReprocessing.unsorted_read_group_quality_distribution_metrics,
                                    ExternalExomeReprocessing.unsorted_read_group_quality_by_cycle_metrics,
                                    ExternalExomeReprocessing.unsorted_read_group_insert_size_metrics,
                                    ExternalExomeReprocessing.unsorted_read_group_base_distribution_by_cycle_metrics,
                                    ExternalExomeReprocessing.quality_yield_metrics,
                                    # File? outputs
                                    select_all([ExternalExomeReprocessing.fingerprint_detail_metrics]),
                                    select_all([ExternalExomeReprocessing.fingerprint_summary_metrics]),
                                    select_all([ExternalExomeReprocessing.cross_check_fingerprints_metrics]),
                                    
    ])

    # Copy results of pipeline to test results bucket
    call Copy.CopyFilesFromCloudToCloud as CopyToTestResults {
      input:
        files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
        vault_token_path          = vault_token_path,
        google_account_vault_path = google_account_vault_path,
        destination_cloud_path    = results_path
    }
  
    # If updating truth then copy output to truth bucket
    if (update_truth){
      call Copy.CopyFilesFromCloudToCloud as CopyToTruth {
        input: 
          files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
          vault_token_path          = vault_token_path,
          google_account_vault_path = google_account_vault_path,
          destination_cloud_path    = truth_path
      }
    }

    # External pipelines only validate by confirming that the reprocessing pipeline succesfully ran
    # and copied the output to the desired destination_cloud_path
    #
    # If the CopyToTestResults step is done then we know that the workflow completed successfully 
    if (!update_truth){

      call VerifyExternalReprocessing.VerifyExternalReprocessing as Verify {
        input:
          done = CopyToTestResults.done,
          destination_cloud_path = destination_cloud_path
      }
    }
}