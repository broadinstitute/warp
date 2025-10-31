version 1.0

# The Arrays pipeline is now deprecated 2025-03-06

import "../../../../pipelines/wdl/genotyping/illumina/IlluminaGenotypingArray.wdl" as IlluminaGenotyping
import "../../../../tasks/wdl/InternalArraysTasks.wdl" as InternalArraysTasks
import "../../../../tasks/wdl/InternalTasks.wdl" as InternalTasks
import "../../../../tasks/wdl/Utilities.wdl" as utils

## Copyright Broad Institute, 2019
##
## This WDL pipeline implements data processing for Illumina Genotyping Arrays
##
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

workflow Arrays {

  String pipeline_version = "2.6.31"

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
    String environment
    File vault_token_path
  }

  if (false) {
    String? none = "None"
  }

  String service_account_filename = "service-account.json"

  # Authorization block to be sourced by tasks that access cloud SQL database
  # Sets up access to vault, reads authentication information from vault, sets permissions for accessing cloud sql.
  Array[String] authentication_block = [
    "export VAULT_ADDR=https://clotho.broadinstitute.org:8200",
    "declare -r secrets=secret/dsde/gotc/~{environment}/metrics/wdl/secrets",
    "for field in user password jdbc_string",
    "do",
    "    n=0",
    "    until [ $n -ge 5 ]",
    "    do",
    "        vault read -field=$field $secrets > cloudsql.db_$field.txt",
    "        if [ -s cloudsql.db_$field.txt ]",
    "        then break",
    "        else",
    "            n=$[$n+1]",
    "            sleep 60",
    "        fi",
    "    done",
    "done",
    "mv cloudsql.db_jdbc_string.txt cloudsql.db_jdbc.txt",
    "if [ ~{environment} == prod ]; then",
    "  key=secret/dsde/gotc/prod/service-accounts/metrics-cloudsql-prod-service-account.json",
    "else",
    "  key=secret/dsde/gotc/dev/service-accounts/metrics-cloudsql-non-prod-service-account.json",
    "fi",
    "vault read --format=json $key | jq .data > ~{service_account_filename}"]

  if ((defined(bead_pool_manifest_filename) || defined(cluster_filename) || defined(gender_cluster_filename) || defined(zcall_thresholds_filename)) && !defined(arrays_metadata_path)) {
    call utils.ErrorWithMessage as ErrorMessageNoArraysMetadataPath {
      input:
        message = "If either bead_pool_manifest_filename, cluster_filename, gender_cluster_filename, or zcall_thresholds_filename is defined, then arrays_metadata_path must also be defined"
    }
  }

  if (!defined(bead_pool_manifest_filename) && !defined(bead_pool_manifest_file)) {
    call utils.ErrorWithMessage as ErrorMessageBpmNoInput {
      input:
        message = "Either bead_pool_manifest_filename or bead_pool_manifest_file (and NOT both) must be defined as input"
    }
  }

  if (defined(bead_pool_manifest_filename) && defined(bead_pool_manifest_file)) {
    call utils.ErrorWithMessage as ErrorMessageBpmDoubleInput {
      input:
        message = "bead_pool_manifest_filename and bead_pool_manifest_file cannot both be defined as input"
    }
  }

  String bpm_filename = if (defined(bead_pool_manifest_file)) then basename(select_first([bead_pool_manifest_file])) else select_first([bead_pool_manifest_filename, ""])
  String chip_type = sub(bpm_filename, ".bpm", "")
  String arrays_chip_metadata_path = select_first([arrays_metadata_path, ""]) + chip_type + "/"
  File bpm_file = if (defined(bead_pool_manifest_file)) then select_first([bead_pool_manifest_file]) else arrays_chip_metadata_path + select_first([bead_pool_manifest_filename, ""])

  if (!defined(cluster_filename) && !defined(cluster_file)) {
    call utils.ErrorWithMessage as ErrorMessageEgtNoInput {
      input:
        message = "Either cluster_filename or cluster_file (and NOT both) must be defined as input"
    }
  }

  if (defined(cluster_filename) && defined(cluster_file)) {
    call utils.ErrorWithMessage as ErrorMessageEgtDoubleInput {
      input:
        message = "cluster_filename and cluster_file cannot both be defined as input"
    }
  }

  File egt_file = if (defined(cluster_file)) then select_first([cluster_file]) else arrays_chip_metadata_path + select_first([cluster_filename, ""])
  String egt_filename = basename(egt_file)

  if (defined(gender_cluster_filename) && defined(gender_cluster_file)) {
    call utils.ErrorWithMessage as ErrorMessageGenderEgtDoubleInput {
      input:
        message = "gender_cluster_filename and gender_cluster_file cannot both be defined as input"
    }
  }

  File? gender_egt_file = if (defined(gender_cluster_file)) then select_first([gender_cluster_file]) else if (defined(gender_cluster_filename)) then arrays_chip_metadata_path + select_first([gender_cluster_filename, ""]) else none

  if (defined(zcall_thresholds_filename) && defined(zcall_thresholds_file)) {
    call utils.ErrorWithMessage as ErrorMessageGenderZCallDoubleInput {
      input:
        message = "zcall_thresholds_filename and zcall_thresholds_file cannot both be defined as input"
    }
  }

  File? zcall_file = if (defined(zcall_thresholds_file)) then select_first([zcall_thresholds_file]) else if (defined(zcall_thresholds_filename)) then arrays_chip_metadata_path + select_first([zcall_thresholds_filename, ""]) else none

  if (!defined(extended_chip_manifest_file)) {
    if (!defined(arrays_metadata_path)) {
      call utils.ErrorWithMessage as ErrorMessageNoExtChipManifestFileNoArraysMetadataPath {
        input:
          message = "If extended_chip_manifest_file is NOT defined, then arrays_metadata_path MUST be defined"
      }
    }
    if (defined(arrays_metadata_path)) {
      String extended_manifest_map_filename = "extended_manifest_map.txt"
      File extended_manifest_map_file = select_first([arrays_metadata_path, ""]) + extended_manifest_map_filename
      call InternalArraysTasks.ResolveExtendedIlluminaManifestFile {
        input:
          extended_manifest_map_file = extended_manifest_map_file,
          bpm_filename = bpm_filename,
          egt_filename = egt_filename,
          arrays_chip_metadata_path = arrays_chip_metadata_path,
          preemptible_tries = preemptible_tries
      }
    }
  }

  if (defined(arrays_metadata_path)) {
    String minor_allele_frequency_map_filename = "maf_map.txt"
    File minor_allele_frequency_map_file = select_first([arrays_metadata_path, ""]) + minor_allele_frequency_map_filename
    call InternalArraysTasks.ResolveMinorAlleleFrequencyFile {
      input:
        minor_allele_frequency_map_file = minor_allele_frequency_map_file,
        bpm_filename = bpm_filename,
        arrays_chip_metadata_path = arrays_chip_metadata_path,
        preemptible_tries = preemptible_tries
    }
    File? resolved_maf = if (ResolveMinorAlleleFrequencyFile.found == true) then ResolveMinorAlleleFrequencyFile.minor_allele_frequency_file else none
  }

  if ((defined(control_sample_name)) &&
      ((!defined(control_sample_vcf_file)) || (!defined(control_sample_vcf_index_file)) || (!defined(control_sample_intervals_file))) &&
      (!defined(arrays_control_data_path))) {
        call utils.ErrorWithMessage as ErrorMessageNoArraysControlDataPath {
          input:
            message = "If either control_sample_name is defined and control_sample_vcf_file, control_sample_vcf_index_file, control_sample_vcf_index, or control_sample_intervals_file ARE NOT defined, then arrays_control_data_path must also be defined"
        }
  }

  File? control_sample_vcf =       if (defined(control_sample_vcf_file))       then select_first([control_sample_vcf_file])       else if (defined(control_sample_name)) then select_first([arrays_control_data_path, ""]) + select_first([control_sample_name, ""]) + ".vcf.gz" else none
  File? control_sample_vcf_index = if (defined(control_sample_vcf_index_file)) then select_first([control_sample_vcf_index_file]) else if (defined(control_sample_name)) then select_first([arrays_control_data_path, ""]) + select_first([control_sample_name, ""]) + ".vcf.gz.tbi" else none
  File? control_sample_intervals = if (defined(control_sample_intervals_file)) then select_first([control_sample_intervals_file]) else if (defined(control_sample_name)) then select_first([arrays_control_data_path, ""]) + select_first([control_sample_name, ""]) + ".interval_list" else none

  if (!defined(params_file)) {
    # If the params_file is not provided, we will generate it from the (currently optional) bunch of parameters.
    # This is to allow for backwards-compatibility.  When we remove params_file as an (optional) input, this will be
    # made a required step and all the inputs will become non-optional
    call InternalArraysTasks.CreateChipWellBarcodeParamsFile {
      input:
        chip_type_name = chip_type,
        chip_well_barcode = chip_well_barcode,
        collaborator_participant_id = collaborator_participant_id,
        lab_batch = lab_batch,
        participant_id = participant_id,
        product_family = product_family,
        product_name = product_name,
        product_order_id = product_order_id,
        product_part_number = product_part_number,
        product_type = product_type,
        regulatory_designation = select_first([regulatory_designation]),
        research_project_id = select_first([research_project_id]),
        sample_alias = sample_alias,
        gender = reported_gender,
        sample_id = sample_id,
        sample_lsid = sample_lsid,
        preemptible_tries = preemptible_tries
    }
  }

  if (!defined(analysis_version_number)) {
    # If the analysis_version_number is not provided as an input, we will query the ARRAYS_QC table to get the
    # next value.
    call InternalArraysTasks.GetNextArraysQcAnalysisVersionNumber {
      input:
        chip_well_barcode = chip_well_barcode,
        preemptible_tries = preemptible_tries,
        vault_token_path = vault_token_path,
        authentication = authentication_block,
        service_account_filename = service_account_filename
    }
  }

  Int analysis_version = select_first([analysis_version_number, GetNextArraysQcAnalysisVersionNumber.analysis_version_number])

  if (read_fingerprint_from_mercury && (!defined(control_sample_name))) {
    call InternalTasks.MakeSafeFilename {
      input:
        name = sample_alias
    }

    call InternalTasks.DownloadGenotypes {
      input:
        sample_alias = sample_alias,
        sample_lsid = sample_lsid,
        output_vcf_base_name = chip_well_barcode + "." + MakeSafeFilename.output_safe_name + ".reference.fingerprint",
        ignoreSpecificGenotypesLsid = sample_lsid,
        ignoreSpecificGenotypesPlatform = "GENERAL_ARRAY",
        haplotype_database_file = haplotype_database_file,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        environment = environment,
        vault_token_path = vault_token_path,
        preemptible_tries = preemptible_tries
    }
  }

  Boolean fingerprint_downloaded_from_mercury = select_first([DownloadGenotypes.fingerprint_retrieved, false])

  call InternalArraysTasks.UpdateChipWellBarcodeIndex {
    input:
      params_file = select_first([CreateChipWellBarcodeParamsFile.params_file, params_file]),
      disk_size = disk_size,
      preemptible_tries = preemptible_tries,
      vault_token_path = vault_token_path,
      authentication = authentication_block,
      service_account_filename = service_account_filename
  }

  File? maf = if (defined(minor_allele_frequency_file)) then select_first([minor_allele_frequency_file]) else resolved_maf

  call IlluminaGenotyping.IlluminaGenotypingArray as IlluminaGenotypingArray {
    input:
      sample_alias = sample_alias,
      analysis_version_number = analysis_version,
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
      bead_pool_manifest_file = bpm_file,
      extended_chip_manifest_file = select_first([extended_chip_manifest_file, ResolveExtendedIlluminaManifestFile.extended_illumina_manifest_file]),
      cluster_file = egt_file,
      gender_cluster_file = gender_egt_file,
      zcall_thresholds_file = zcall_file,
      fingerprint_genotypes_vcf_file = if (fingerprint_downloaded_from_mercury) then DownloadGenotypes.reference_fingerprint_vcf else fingerprint_genotypes_vcf_file,
      fingerprint_genotypes_vcf_index_file = if (fingerprint_downloaded_from_mercury) then DownloadGenotypes.reference_fingerprint_vcf_index else fingerprint_genotypes_vcf_index_file,
      haplotype_database_file = haplotype_database_file,
      variant_rsids_file = variant_rsids_file,
      subsampled_metrics_interval_list = subsampled_metrics_interval_list,
      contamination_controls_vcf = contamination_controls_vcf,
      minor_allele_frequency_file = maf,
      control_sample_vcf_file = control_sample_vcf,
      control_sample_vcf_index_file = control_sample_vcf_index,
      control_sample_intervals_file = control_sample_intervals,
      control_sample_name = control_sample_name,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries,
      genotype_concordance_threshold = genotype_concordance_threshold
  }

  if (size(IlluminaGenotypingArray.gtc) == 0) {
    call InternalArraysTasks.GenerateEmptyVariantCallingMetricsFile {
      input:
        chip_well_barcode = chip_well_barcode,
        sample_alias = sample_alias,
        chip_type = chip_type,
        reported_gender = reported_gender,
        autocall_version = IlluminaGenotypingArray.autocall_version,
        output_metrics_basename = sample_alias,
        cluster_filename = egt_filename,
        analysis_version_number = analysis_version,
        preemptible_tries = preemptible_tries
    }

    call InternalArraysTasks.UploadEmptyArraysMetrics {
      input:
        arrays_variant_calling_detail_metrics = GenerateEmptyVariantCallingMetricsFile.detail_metrics,
        disk_size = disk_size,
        preemptible_tries = preemptible_tries,
        vault_token_path = vault_token_path,
        authentication = authentication_block,
        service_account_filename = service_account_filename
    }

    call InternalArraysTasks.BlacklistBarcode as BlacklistFailedNormalization {
      input:
        upload_metrics_output = UploadEmptyArraysMetrics.upload_metrics_empty_file,
        analysis_version_number = analysis_version,
        chip_well_barcode = chip_well_barcode,
        preemptible_tries = preemptible_tries,
        reason = "DATA_QUALITY",
        notes = "Normalization Failed",
        vault_token_path = vault_token_path,
        authentication = authentication_block,
        service_account_filename = service_account_filename
    }
  }

  if (size(IlluminaGenotypingArray.gtc) > 0) {
    call InternalArraysTasks.VcfToMercuryFingerprintJson {
      input:
        input_vcf_file = select_first([IlluminaGenotypingArray.output_fingerprint_vcf]),
        input_vcf_index_file = select_first([IlluminaGenotypingArray.output_fingerprint_vcf_index]),
        variant_calling_detail_metrics_file = select_first([IlluminaGenotypingArray.arrays_variant_calling_detail_metrics]),
        sample_lsid = sample_lsid,
        output_json_filename = chip_well_barcode + ".fingerprint.json",
        disk_size = disk_size,
        preemptible_tries = preemptible_tries
    }

    # Only write fingerprints to the Mercury Fingerprint Store if writing is enabled AND the sample is NOT a control
    if (write_fingerprint_to_mercury && (!defined(control_sample_name))) {
      call InternalTasks.UploadFingerprintToMercury {
        input:
          fingerprint_json_file = VcfToMercuryFingerprintJson.output_json_file,
          gtc_file = IlluminaGenotypingArray.gtc,
          environment = environment,
          vault_token_path = vault_token_path,
          preemptible_tries = preemptible_tries
      }
    }

    if (defined(IlluminaGenotypingArray.bafregress_results_file)) {
      call InternalArraysTasks.CreateBafRegressMetricsFile {
        input:
          input_file = select_first([IlluminaGenotypingArray.bafregress_results_file]),
          output_metrics_basefilename = chip_well_barcode,
          disk_size = disk_size,
          preemptible_tries = preemptible_tries
      }
    }

    call InternalArraysTasks.UploadArraysMetrics {
      input:
        arrays_variant_calling_detail_metrics = select_first([IlluminaGenotypingArray.arrays_variant_calling_detail_metrics]),
        arrays_variant_calling_summary_metrics = select_first([IlluminaGenotypingArray.arrays_variant_calling_summary_metrics]),
        arrays_control_code_summary_metrics = select_first([IlluminaGenotypingArray.arrays_variant_calling_control_metrics]),
        fingerprinting_detail_metrics = IlluminaGenotypingArray.fingerprint_detail_metrics,
        fingerprinting_summary_metrics = IlluminaGenotypingArray.fingerprint_summary_metrics,
        genotype_concordance_summary_metrics = IlluminaGenotypingArray.genotype_concordance_summary_metrics,
        genotype_concordance_detail_metrics  = IlluminaGenotypingArray.genotype_concordance_detail_metrics,
        genotype_concordance_contingency_metrics = IlluminaGenotypingArray.genotype_concordance_contingency_metrics,
        verify_id_metrics = IlluminaGenotypingArray.contamination_metrics,
        bafregress_metrics = CreateBafRegressMetricsFile.output_metrics_file,
        disk_size = disk_size,
        preemptible_tries = preemptible_tries,
        vault_token_path = vault_token_path,
        authentication = authentication_block,
        service_account_filename = service_account_filename
    }

    if (select_first([IlluminaGenotypingArray.genotype_concordance_failed, false])) {
      call InternalArraysTasks.BlacklistBarcode as BlacklistFailedGenotypeConcordance {
        input:
          upload_metrics_output = UploadArraysMetrics.upload_metrics_empty_file,
          analysis_version_number = analysis_version,
          chip_well_barcode = chip_well_barcode,
          preemptible_tries = preemptible_tries,
          reason = "GENOTYPE_CONCORDANCE",
          notes = "Genotype concordance below threshold: ~{genotype_concordance_threshold}",
          vault_token_path = vault_token_path,
          authentication = authentication_block,
          service_account_filename = service_account_filename
      }
    }

    if (defined(fingerprint_genotypes_vcf_file) &&
        select_first([IlluminaGenotypingArray.check_fingerprint_lod, 0.0]) < -3.0) {
      call InternalArraysTasks.BlacklistBarcode as BlacklistFailedFingerprint {
        input:
          upload_metrics_output = UploadArraysMetrics.upload_metrics_empty_file,
          analysis_version_number = analysis_version,
          chip_well_barcode = chip_well_barcode,
          preemptible_tries = preemptible_tries,
          reason = "SAMPLE_MIXUP",
          notes = "Fingerprint LOD below -3.0",
          vault_token_path = vault_token_path,
          authentication = authentication_block,
          service_account_filename = service_account_filename
      }
    }
  }

  output {
    String chip_well_barcode_output = IlluminaGenotypingArray.chip_well_barcode_output
    Int analysis_version_number_output = IlluminaGenotypingArray.analysis_version_number_output
    File gtc_file = IlluminaGenotypingArray.gtc
    File red_idat_md5_cloud_path = IlluminaGenotypingArray.red_idat_md5_cloud_path
    File green_idat_md5_cloud_path = IlluminaGenotypingArray.green_idat_md5_cloud_path
    File output_bead_pool_manifest_file = bpm_file
    File? output_vcf_md5_cloud_path = IlluminaGenotypingArray.output_vcf_md5_cloud_path
    File? output_vcf = IlluminaGenotypingArray.output_vcf
    File? output_vcf_index = IlluminaGenotypingArray.output_vcf_index
    File? baf_regress_metrics_file = CreateBafRegressMetricsFile.output_metrics_file
    File? contamination_metrics_file = IlluminaGenotypingArray.contamination_metrics
    File? reference_fingerprint_vcf = DownloadGenotypes.reference_fingerprint_vcf
    File? reference_fingerprint_vcf_index = DownloadGenotypes.reference_fingerprint_vcf_index
    File? output_fingerprint_vcf = IlluminaGenotypingArray.output_fingerprint_vcf
    File? output_fingerprint_vcf_index = IlluminaGenotypingArray.output_fingerprint_vcf_index
    File? output_fingerprint_json_file = VcfToMercuryFingerprintJson.output_json_file
    File arrays_variant_calling_detail_metrics_file = select_first([IlluminaGenotypingArray.arrays_variant_calling_detail_metrics, GenerateEmptyVariantCallingMetricsFile.detail_metrics])
    File? arrays_variant_calling_summary_metrics_file = IlluminaGenotypingArray.arrays_variant_calling_summary_metrics
    File? arrays_variant_calling_control_metrics_file = IlluminaGenotypingArray.arrays_variant_calling_control_metrics
    File? arrays_subset_variant_calling_detail_metrics_file = IlluminaGenotypingArray.arrays_subset_variant_calling_detail_metrics
    File? arrays_subset_variant_calling_summary_metrics_file = IlluminaGenotypingArray.arrays_subset_variant_calling_summary_metrics
    File? arrays_subset_variant_calling_control_metrics_file = IlluminaGenotypingArray.arrays_subset_variant_calling_control_metrics
    File? fingerprint_detail_metrics_file = IlluminaGenotypingArray.fingerprint_detail_metrics
    File? fingerprint_summary_metrics_file = IlluminaGenotypingArray.fingerprint_summary_metrics
    File? genotype_concordance_summary_metrics_file = IlluminaGenotypingArray.genotype_concordance_summary_metrics
    File? genotype_concordance_detail_metrics_file  = IlluminaGenotypingArray.genotype_concordance_detail_metrics
    File? genotype_concordance_contingency_metrics_file = IlluminaGenotypingArray.genotype_concordance_contingency_metrics
    File chip_well_barcode_params_file = select_first([CreateChipWellBarcodeParamsFile.params_file, params_file])
  }
  meta {
    allowNestedInputs: true
  }
}
