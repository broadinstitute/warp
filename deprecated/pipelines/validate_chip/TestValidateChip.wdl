version 1.0


import "../../pipelines/wdl/arrays/validate_chip/ValidateChip.wdl" as ValidateChip
import "../../verification/VerifyValidateChip.wdl" as VerifyValidateChip
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/CopyFilesFromCloudToCloud.wdl" as Copy

workflow TestValidateChip {

    input {
      String sample_alias
      Int analysis_version_number
      Float call_rate_threshold
      String reported_gender
      String chip_well_barcode
      File red_idat_cloud_path
      File green_idat_cloud_path
      File ref_fasta
      File ref_fasta_index
      File ref_dict
      File dbSNP_vcf
      File dbSNP_vcf_index
      File bead_pool_manifest_file
      String chip_type = basename(bead_pool_manifest_file, ".bpm")
      File chip_manifest_csv_file
      File supported_ref_fasta
      File supported_ref_fasta_index
      File supported_ref_dict
      File chain_file
      File cluster_file
      File control_sample_vcf_file
      File control_sample_vcf_index_file
      File control_sample_intervals_file
      String control_sample_name
      Float indel_genotype_concordance_threshold
      Int disk_size
      Int preemptible_tries

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
  
    call ValidateChip.ValidateChip {
      input:
        sample_alias = sample_alias,
        analysis_version_number = analysis_version_number,
        call_rate_threshold = call_rate_threshold,
        reported_gender = reported_gender,
        chip_well_barcode = chip_well_barcode,
        red_idat_cloud_path = red_idat_cloud_path,
        green_idat_cloud_path = green_idat_cloud_path,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        bead_pool_manifest_file = bead_pool_manifest_file,
        chip_type = chip_type,
        chip_manifest_csv_file = chip_manifest_csv_file,
        supported_ref_fasta = supported_ref_fasta,
        supported_ref_fasta_index = supported_ref_fasta_index,
        supported_ref_dict = supported_ref_dict,
        chain_file = chain_file,
        cluster_file = cluster_file,
        control_sample_vcf_file = control_sample_vcf_file,
        control_sample_vcf_index_file = control_sample_vcf_index_file,
        control_sample_intervals_file = control_sample_intervals_file,
        control_sample_name = control_sample_name,
        indel_genotype_concordance_threshold = indel_genotype_concordance_threshold,
        disk_size = disk_size,
        preemptible_tries = preemptible_tries
  
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    ValidateChip.indel_genotype_concordance_txt_file,
                                    ValidateChip.indel_genotype_concordance_vcf,
                                    ValidateChip.genotype_concordance_txt_file,
                                    ValidateChip.genotype_concordance_vcf,
                                    ValidateChip.create_extended_illumina_manifest_report_file,
                                    ValidateChip.create_extended_illumina_manifest_bad_assays_file,
                                    ValidateChip.create_extended_illumina_manifest_extended_csv_file,
                                    ValidateChip.output_vcf_index,
                                    ValidateChip.output_vcf,
                                    ValidateChip.gtc_file,
                                    ],
                                    
    ])

    
    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    [ # File outputs
                                    ValidateChip.indel_genotype_concordance_contingency_metrics_file,
                                    ValidateChip.indel_genotype_concordance_detail_metrics_file,
                                    ValidateChip.indel_genotype_concordance_summary_metrics_file,
                                    ValidateChip.genotype_concordance_contingency_metrics_file,
                                    ValidateChip.genotype_concordance_detail_metrics_file,
                                    ValidateChip.genotype_concordance_summary_metrics_file,
                                    ValidateChip.arrays_variant_calling_control_metrics_file,
                                    ValidateChip.arrays_variant_calling_summary_metrics_file,
                                    ValidateChip.arrays_variant_calling_detail_metrics_file,
                                    ],
                                    
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

    # If not updating truth then we need to collect all input for the validation WDL
    # This is achieved by passing each desired file/array[files] to GetValidationInputs
    if (!update_truth){
        call Utilities.GetValidationInputs as GetGtc {
          input:
            input_file = ValidateChip.gtc_file,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetVcf {
          input:
            input_file = ValidateChip.output_vcf,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetGenotypeConcordanceVcf {
          input:
            input_file = ValidateChip.genotype_concordance_vcf,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetIndelGenotypeConcordanceVcf {
          input:
            input_file = ValidateChip.indel_genotype_concordance_vcf,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetMetrics {
          input:
            input_files = pipeline_metrics,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifyValidateChip.VerifyValidateChip as Verify {
        input:
          truth_gtc = GetGtc.truth_file, 
          test_gtc = GetGtc.results_file,
          truth_vcf = GetVcf.truth_file, 
          test_vcf = GetVcf.results_file,
          truth_genotype_concordance_vcf = GetGenotypeConcordanceVcf.truth_file, 
          test_genotype_concordance_vcf = GetGenotypeConcordanceVcf.results_file,
          truth_indel_genotype_concordance_vcf = GetIndelGenotypeConcordanceVcf.truth_file, 
          test_indel_genotype_concordance_vcf = GetIndelGenotypeConcordanceVcf.results_file,
          truth_metrics = GetMetrics.truth_files, 
          test_metrics = GetMetrics.results_files,
          done = CopyToTestResults.done,
          bead_pool_manifest_file = ValidateChip.output_bead_pool_manifest_file
      }

    }

    output {}

}