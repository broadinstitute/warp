version 1.0


import "../../pipelines/wdl/reprocessing/external/wgs/ExternalWholeGenomeReprocessing.wdl" as ExternalWholeGenomeReprocessing
import "../../verification/VerifyExternalReprocessing.wdl" as VerifyExternalReprocessing
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/CopyFilesFromCloudToCloud.wdl" as Copy

workflow TestExternalWholeGenomeReprocessing {

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
  
    call ExternalWholeGenomeReprocessing.ExternalWholeGenomeReprocessing {
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
        wgs_coverage_interval_list = wgs_coverage_interval_list,
        fingerprint_genotypes_file = fingerprint_genotypes_file,
        fingerprint_genotypes_index = fingerprint_genotypes_index,
        references = references,
        scatter_settings = scatter_settings,
        papi_settings = papi_settings,
        destination_cloud_path = destination_cloud_path,
        vault_token_path = vault_token_path,
        google_account_vault_path = google_account_vault_path,
        cloud_provider = cloud_provider
  
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    ExternalWholeGenomeReprocessing.output_vcf_index,
                                    ExternalWholeGenomeReprocessing.output_vcf,
                                    ExternalWholeGenomeReprocessing.validate_cram_file_report,
                                    ExternalWholeGenomeReprocessing.output_cram_md5,
                                    ExternalWholeGenomeReprocessing.output_cram_index,
                                    ExternalWholeGenomeReprocessing.output_cram,
                                    ExternalWholeGenomeReprocessing.agg_quality_distribution_pdf,
                                    ExternalWholeGenomeReprocessing.agg_insert_size_histogram_pdf,
                                    ExternalWholeGenomeReprocessing.agg_gc_bias_pdf,
                                    ExternalWholeGenomeReprocessing.calculate_read_group_checksum_md5,
                                    ExternalWholeGenomeReprocessing.selfSM,
                                    ExternalWholeGenomeReprocessing.read_group_gc_bias_pdf,
                                    ],
                                    # Array[File] outputs
                                    ExternalWholeGenomeReprocessing.unsorted_read_group_quality_distribution_pdf,
                                    ExternalWholeGenomeReprocessing.unsorted_read_group_quality_by_cycle_pdf,
                                    ExternalWholeGenomeReprocessing.unsorted_read_group_insert_size_histogram_pdf,
                                    ExternalWholeGenomeReprocessing.unsorted_read_group_base_distribution_by_cycle_pdf,
                                    ExternalWholeGenomeReprocessing.unmapped_bams,
                                    ExternalWholeGenomeReprocessing.validation_report,
                                    # File? outputs
                                    select_all([ExternalWholeGenomeReprocessing.output_bqsr_reports]),
                                    
    ])

    
    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    [ # File outputs
                                    ExternalWholeGenomeReprocessing.gvcf_detail_metrics,
                                    ExternalWholeGenomeReprocessing.gvcf_summary_metrics,
                                    ExternalWholeGenomeReprocessing.duplicate_metrics,
                                    ExternalWholeGenomeReprocessing.raw_wgs_metrics,
                                    ExternalWholeGenomeReprocessing.wgs_metrics,
                                    ExternalWholeGenomeReprocessing.agg_quality_distribution_metrics,
                                    ExternalWholeGenomeReprocessing.agg_pre_adapter_summary_metrics,
                                    ExternalWholeGenomeReprocessing.agg_pre_adapter_detail_metrics,
                                    ExternalWholeGenomeReprocessing.agg_insert_size_metrics,
                                    ExternalWholeGenomeReprocessing.agg_gc_bias_summary_metrics,
                                    ExternalWholeGenomeReprocessing.agg_gc_bias_detail_metrics,
                                    ExternalWholeGenomeReprocessing.agg_bait_bias_summary_metrics,
                                    ExternalWholeGenomeReprocessing.agg_bait_bias_detail_metrics,
                                    ExternalWholeGenomeReprocessing.agg_alignment_summary_metrics,
                                    ExternalWholeGenomeReprocessing.read_group_gc_bias_summary_metrics,
                                    ExternalWholeGenomeReprocessing.read_group_gc_bias_detail_metrics,
                                    ExternalWholeGenomeReprocessing.read_group_alignment_summary_metrics,
                                    ],
                                    # Array[File] outputs
                                    ExternalWholeGenomeReprocessing.unsorted_read_group_quality_distribution_metrics,
                                    ExternalWholeGenomeReprocessing.unsorted_read_group_quality_by_cycle_metrics,
                                    ExternalWholeGenomeReprocessing.unsorted_read_group_insert_size_metrics,
                                    ExternalWholeGenomeReprocessing.unsorted_read_group_base_distribution_by_cycle_metrics,
                                    ExternalWholeGenomeReprocessing.quality_yield_metrics,
                                    # File? outputs
                                    select_all([ExternalWholeGenomeReprocessing.fingerprint_detail_metrics]),
                                    select_all([ExternalWholeGenomeReprocessing.fingerprint_summary_metrics]),
                                    select_all([ExternalWholeGenomeReprocessing.cross_check_fingerprints_metrics]),
                                    
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

    # This is achieved by passing each desired file/array[files] to GetValidationInputs
    if (!update_truth){

      call VerifyExternalReprocessing.VerifyExternalReprocessing as Verify {
        input:
          done = CopyToTestResults.done,
          destination_cloud_path = destination_cloud_path
      }
    }
}