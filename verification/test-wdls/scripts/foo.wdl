version 1.0


import "../../pipelines/broad/genotyping/illumina/IlluminaGenotypingArray.wdl" as IlluminaGenotyping
import "../../tasks/broad/InternalArraysTasks.wdl" as InternalArraysTasks
import "../../tasks/broad/InternalTasks.wdl" as InternalTasks
import "../../tasks/broad/Utilities.wdl" as utils
import "../../pipelines/broad/arrays/single_sample/Arrays.wdl" as Arrays

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
      Boolean read_fingerprint_from_mercury = false
      File? fingerprint_genotypes_vcf_file
      File? fingerprint_genotypes_vcf_index_file
      File haplotype_database_file
      File variant_rsids_file
      Boolean write_fingerprint_to_mercury = false
      File? subsampled_metrics_interval_list
      File? contamination_controls_vcf
      File? minor_allele_frequency_file
      String? arrays_control_data_path
      String? control_sample_name
      File? control_sample_vcf_file
      File? control_sample_vcf_index_file
      File? control_sample_intervals_file
      Int disk_size
      Int preemptible_tries
      String environment
      File vault_token_path

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

  call Arrays.Arrays {
    input:
      chip_well_barcode = chip_well_barcode,
      analysis_version_number = analysis_version_number,
      sample_alias = sample_alias,
      sample_id = sample_id,
      sample_lsid = sample_lsid,
      reported_gender = reported_gender,
      collaborator_participant_id = collaborator_participant_id,
      participant_id = participant_id,
      call_rate_threshold = call_rate_threshold,
      genotype_concordance_threshold = genotype_concordance_threshold,
      red_idat_cloud_path = red_idat_cloud_path,
      green_idat_cloud_path = green_idat_cloud_path,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      params_file = params_file,
      lab_batch = lab_batch,
      product_family = product_family,
      product_name = product_name,
      product_order_id = product_order_id,
      product_part_number = product_part_number,
      product_type = product_type,
      regulatory_designation = regulatory_designation,
      research_project_id = research_project_id,
      arrays_metadata_path = arrays_metadata_path,
      bead_pool_manifest_filename = bead_pool_manifest_filename,
      bead_pool_manifest_file = bead_pool_manifest_file,
      cluster_filename = cluster_filename,
      cluster_file = cluster_file,
      gender_cluster_filename = gender_cluster_filename,
      gender_cluster_file = gender_cluster_file,
      zcall_thresholds_filename = zcall_thresholds_filename,
      zcall_thresholds_file = zcall_thresholds_file,
      extended_chip_manifest_file = extended_chip_manifest_file,
      read_fingerprint_from_mercury = read_fingerprint_from_mercury,
      fingerprint_genotypes_vcf_file = fingerprint_genotypes_vcf_file,
      fingerprint_genotypes_vcf_index_file = fingerprint_genotypes_vcf_index_file,
      haplotype_database_file = haplotype_database_file,
      variant_rsids_file = variant_rsids_file,
      write_fingerprint_to_mercury = write_fingerprint_to_mercury,
      subsampled_metrics_interval_list = subsampled_metrics_interval_list,
      contamination_controls_vcf = contamination_controls_vcf,
      minor_allele_frequency_file = minor_allele_frequency_file,
      arrays_control_data_path = arrays_control_data_path,
      control_sample_name = control_sample_name,
      control_sample_vcf_file = control_sample_vcf_file,
      control_sample_vcf_index_file = control_sample_vcf_index_file,
      control_sample_intervals_file = control_sample_intervals_file,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries,
      environment = environment,
      vault_token_path = vault_token_path

  }

  # Collect all of the pipeline outputs into single Array[String]
  Array[String] pipeline_outputs = flatten([

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



}