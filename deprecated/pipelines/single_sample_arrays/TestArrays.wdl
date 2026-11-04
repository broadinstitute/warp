version 1.0

import "../../pipelines/wdl/arrays/single_sample/Arrays.wdl" as Arrays
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/CopyFilesFromCloudToCloud.wdl" as Copy
import "../../verification/VerifyArrays.wdl" as VerifyArrays

workflow TestArrays {

    input {
      String chip_well_barcode
      Int? analysis_version_number
  
      String sample_alias
      String? sample_id
      String sample_lsid
      String reported_gender
      String? collaborator_participant_id
      String? participant_id
  
      Float call_rate_threshold = 0.98
      Float genotype_concordance_threshold = 0.95
  
      File red_idat_cloud_path
      File green_idat_cloud_path
      File ref_fasta
      File ref_fasta_index
      File ref_dict
  
      File dbSNP_vcf
      File dbSNP_vcf_index
  
      File? params_file
      String? lab_batch
      String? product_family
      String? product_name
      String? product_order_id
      String? product_part_number
      String product_type = ""
      String? regulatory_designation
      String? research_project_id
  
      String? arrays_metadata_path
  
      String? bead_pool_manifest_filename
      File? bead_pool_manifest_file
  
      String? cluster_filename
      File? cluster_file
  
      String? gender_cluster_filename
      File? gender_cluster_file
  
      String? zcall_thresholds_filename
      File? zcall_thresholds_file
  
      File? extended_chip_manifest_file
  
      # For CheckFingerprint:
      # If this is true, we will read fingerprints from Mercury
      # Otherwise, we will use the optional input fingerprint VCFs below
      Boolean read_fingerprint_from_mercury = false
      File? fingerprint_genotypes_vcf_file
      File? fingerprint_genotypes_vcf_index_file
      File haplotype_database_file
  
      # For fingerprint generation.
      # For SelectVariants, in order to generate a fingerprint for upload to Mercury
      File variant_rsids_file
      # If this is true, the WDL will upload the generated finerprint to Mercury
      # The generated fingerprint VCF and json are available as outputs too.
      Boolean write_fingerprint_to_mercury = false
  
      # For Subsampled Metrics
      File? subsampled_metrics_interval_list
  
      # For Contamination Checking
      File? contamination_controls_vcf
  
      # For BAFRegress
      File? minor_allele_frequency_file
  
      # For HapMap GenotypeConcordance Check:
      String? arrays_control_data_path
  
      String? control_sample_name
      File? control_sample_vcf_file
      File? control_sample_vcf_index_file
      File? control_sample_intervals_file
  
      Int disk_size
      Int preemptible_tries
    
      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
      String vault_token_path
      String vault_token_path_arrays
      String google_account_vault_path
      String environment

    }

    call Arrays.Arrays {
      input:
        chip_well_barcode                    = chip_well_barcode,
        analysis_version_number              = analysis_version_number,
        sample_alias                         = sample_alias,
        sample_id                            = sample_id,
        sample_lsid                          = sample_lsid,
        reported_gender                      = reported_gender,
        collaborator_participant_id          = collaborator_participant_id,
        participant_id                       = participant_id,
        call_rate_threshold                  = call_rate_threshold,
        genotype_concordance_threshold       = genotype_concordance_threshold,
        red_idat_cloud_path                  = red_idat_cloud_path,
        green_idat_cloud_path                = green_idat_cloud_path,
        ref_fasta                            = ref_fasta,
        ref_fasta_index                      = ref_fasta_index,
        ref_dict                             = ref_dict,
        dbSNP_vcf                            = dbSNP_vcf,
        dbSNP_vcf_index                      = dbSNP_vcf_index,
        params_file                          = params_file,
        lab_batch                            = lab_batch,
        product_family                       = product_family,
        product_name                         = product_name,
        product_order_id                     = product_order_id,
        product_part_number                  = product_part_number,
        product_type                         = product_type,
        regulatory_designation               = regulatory_designation,
        research_project_id                  = research_project_id,
        arrays_metadata_path                 = arrays_metadata_path,
        bead_pool_manifest_filename          = bead_pool_manifest_filename,
        bead_pool_manifest_file              = bead_pool_manifest_file,
        cluster_filename                     = cluster_filename,
        cluster_file                         = cluster_file,
        gender_cluster_filename              = gender_cluster_filename,
        gender_cluster_file                  = gender_cluster_file,
        zcall_thresholds_filename            = zcall_thresholds_filename,
        zcall_thresholds_file                = zcall_thresholds_file,
        extended_chip_manifest_file          = extended_chip_manifest_file,
        read_fingerprint_from_mercury        = read_fingerprint_from_mercury,
        fingerprint_genotypes_vcf_file       = fingerprint_genotypes_vcf_file,
        fingerprint_genotypes_vcf_index_file = fingerprint_genotypes_vcf_index_file,
        haplotype_database_file              = haplotype_database_file,
        variant_rsids_file                   = variant_rsids_file,
        write_fingerprint_to_mercury         = write_fingerprint_to_mercury,
        subsampled_metrics_interval_list     = subsampled_metrics_interval_list,
        contamination_controls_vcf           = contamination_controls_vcf,
        minor_allele_frequency_file          = minor_allele_frequency_file,
        arrays_control_data_path             = arrays_control_data_path,
        control_sample_name                  = control_sample_name,
        control_sample_vcf_file              = control_sample_vcf_file,
        control_sample_vcf_index_file        = control_sample_vcf_index_file,
        control_sample_intervals_file        = control_sample_intervals_file,
        disk_size                            = disk_size,
        preemptible_tries                    = preemptible_tries,
        environment                          = environment,
        vault_token_path                     = vault_token_path_arrays
    }

    # Collect all of the pipeline outputs into singe Array[String]
    Array[String] pipeline_outputs = flatten([
                            [ # File outputs
                            Arrays.gtc_file,
                            Arrays.red_idat_md5_cloud_path,
                            Arrays.green_idat_md5_cloud_path,
                            Arrays.chip_well_barcode_params_file,
                            Arrays.output_bead_pool_manifest_file
                            ], # File? outputs
                            select_all([Arrays.output_vcf_md5_cloud_path]),
                            select_all([Arrays.output_vcf]),
                            select_all([Arrays.output_vcf_index]),
                            select_all([Arrays.reference_fingerprint_vcf]),
                            select_all([Arrays.reference_fingerprint_vcf_index]),
                            select_all([Arrays.output_fingerprint_vcf]),
                            select_all([Arrays.output_fingerprint_vcf_index]),
                            select_all([Arrays.output_fingerprint_json_file]),
    ])

    # Collect all of the pipeline metric into a single Array[String]
    Array[String] pipeline_metrics = flatten([
                            [ # File_outputs
                            Arrays.arrays_variant_calling_detail_metrics_file
                            ], # File? outputs
                            select_all([Arrays.baf_regress_metrics_file]),
                            select_all([Arrays.contamination_metrics_file]),
                            select_all([Arrays.arrays_variant_calling_summary_metrics_file]),
                            select_all([Arrays.arrays_variant_calling_control_metrics_file]),
                            select_all([Arrays.arrays_subset_variant_calling_detail_metrics_file]),
                            select_all([Arrays.arrays_subset_variant_calling_summary_metrics_file]),
                            select_all([Arrays.arrays_subset_variant_calling_control_metrics_file]),
                            select_all([Arrays.fingerprint_detail_metrics_file]),
                            select_all([Arrays.fingerprint_summary_metrics_file]),
                            select_all([Arrays.genotype_concordance_summary_metrics_file]),
                            select_all([Arrays.genotype_concordance_detail_metrics_file]),
                            select_all([Arrays.genotype_concordance_contingency_metrics_file]),
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
          destination_cloud_path    = results_path
      }
    }

    # If not updating truth then we need to collect all input for the validation WDL
    # This is achieved by passing each desired file/array[files] to GetValidationInputs
    if (!update_truth) {
      call Utilities.GetValidationInputs as GetMetricsInputs {
        input:
          input_files  = pipeline_metrics,
          results_path = results_path,
          truth_path   = truth_path
        }

      call Utilities.GetValidationInputs as GetGtcs {
        input:
            input_file   = Arrays.gtc_file,
            results_path = results_path,
            truth_path   = truth_path
        }
  
      call Utilities.GetValidationInputs as GetVCFs {
        input:
          input_file   = Arrays.output_vcf,
          results_path = results_path,
          truth_path   = truth_path
      }
  
      call Utilities.GetValidationInputs as GetFingerprintVCFs { 
        input:
          input_file   = Arrays.output_fingerprint_vcf,
          results_path = results_path,
          truth_path   = truth_path
      }

      call Utilities.GetValidationInputs as GetRedIdat { 
        input:
          input_file   = Arrays.red_idat_md5_cloud_path,
          results_path = results_path,
          truth_path   = truth_path
      }

      call Utilities.GetValidationInputs as GetGreenIdat { 
        input:
          input_file   = Arrays.green_idat_md5_cloud_path,
          results_path = results_path,
          truth_path   = truth_path
      }

      call Utilities.GetValidationInputs as GetParams { 
        input:
          input_file   = Arrays.chip_well_barcode_params_file,
          results_path = results_path,
          truth_path   = truth_path
      }

      call VerifyArrays.VerifyArrays as VerifyArrays {
        input:
          truth_metrics           = GetMetricsInputs.truth_files,
          truth_gtc               = GetGtcs.truth_file,
          truth_vcf               = GetVCFs.truth_file,
          truth_fp_vcf            = GetFingerprintVCFs.truth_file,
          truth_red_idat_md5      = GetRedIdat.truth_file,
          truth_green_idat_md5    = GetGreenIdat.truth_file,
          truth_params_file       = GetParams.truth_file,
          test_metrics            = GetMetricsInputs.results_files,
          test_gtc                = GetGtcs.results_file,
          test_vcf                = GetVCFs.results_file,
          test_fp_vcf             = GetFingerprintVCFs.results_file,
          test_red_idat_md5       = GetRedIdat.results_file,
          test_green_idat_md5     = GetGreenIdat.results_file,
          test_params_file        = GetParams.results_file,
          bead_pool_manifest_file = Arrays.output_bead_pool_manifest_file,
          done                    = CopyToTestResults.done
        }

    }


}