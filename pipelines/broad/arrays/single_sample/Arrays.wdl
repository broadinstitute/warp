version 1.0

import "../../../../pipelines/broad/genotyping/illumina/IlluminaGenotypingArray.wdl" as IlluminaGenotyping
import "../../../../tasks/broad/InternalArraysTasks.wdl" as InternalTasks

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

  String pipeline_version = "2.3.1"

  input {

    # This is the autocall_version, needed for the case where autocall fails (likely due to normalization errors)
    # In this case it no longer emits the version in its output, so we store it here.
    String autocall_version = "3.0.0"
    String sample_alias
    String sample_lsid
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

    File params_file

    File bead_pool_manifest_file
    String chip_type = basename(bead_pool_manifest_file, ".bpm")

    File extended_chip_manifest_file
    File cluster_file
    String cluster_filename = basename(cluster_file)
    File? gender_cluster_file
    File? zcall_thresholds_file

    # For CheckFingerprint:
    File? fingerprint_genotypes_vcf_file
    File? fingerprint_genotypes_vcf_index_file
    File haplotype_database_file

    # For SelectVariants
    File variant_rsids_file

    # For Subsampled Metrics
    File? subsampled_metrics_interval_list

    # For Contamination Checking
    File? contamination_controls_vcf

    # For BAFRegress
    File? minor_allele_frequency_file

    # For HapMap GenotypeConcordance Check:
    File? control_sample_vcf_file
    File? control_sample_vcf_index_file
    File? control_sample_intervals_file
    String? control_sample_name

    Int disk_size
    Int preemptible_tries
    String environment
    File vault_token_path

    Float genotype_concordance_threshold = 0.98
  }

  String service_account_filename = "service-account.json"

  Array[String] vault_source_block = ["export VAULT_ADDR=https://clotho.broadinstitute.org:8200",
                                      "export VAULT_TOKEN=~{read_lines(vault_token_path)[0]}"]

  # An array of strings to be run by VMs that need authentication
  Array[String] authentication_block = ["export VAULT_ADDR=https://clotho.broadinstitute.org:8200",
  "export VAULT_TOKEN=~{read_lines(vault_token_path)[0]}",
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

  call DownloadGenotypes {
    input:
      sample_alias = sample_alias,
      sample_lsid = sample_lsid,
      output_vcf_name = chip_well_barcode + "." + sample_alias + ".reference.fingerprint.gz",
      haplotype_database_file = haplotype_database_file,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      preemptible_tries = preemptible_tries,
      source_block = vault_source_block
  }

  call InternalTasks.UpdateChipWellBarcodeIndex {
    input:
      params_file = params_file,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries,
      authentication = authentication_block,
      service_account_filename = service_account_filename
  }

  call IlluminaGenotyping.IlluminaGenotypingArray as IlluminaGenotypingArray {
    input:
      autocall_version = autocall_version,
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
      extended_chip_manifest_file = extended_chip_manifest_file,
      cluster_file = cluster_file,
      gender_cluster_file = gender_cluster_file,
      zcall_thresholds_file = zcall_thresholds_file,
      fingerprint_genotypes_vcf_file = DownloadGenotypes.output_vcf,
      fingerprint_genotypes_vcf_index_file = DownloadGenotypes.output_vcf_index,
      haplotype_database_file = haplotype_database_file,
      variant_rsids_file = variant_rsids_file,
      subsampled_metrics_interval_list = subsampled_metrics_interval_list,
      contamination_controls_vcf = contamination_controls_vcf,
      minor_allele_frequency_file = minor_allele_frequency_file,
      control_sample_vcf_file = control_sample_vcf_file,
      control_sample_vcf_index_file = control_sample_vcf_index_file,
      control_sample_intervals_file = control_sample_intervals_file,
      control_sample_name = control_sample_name,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries,
      genotype_concordance_threshold = genotype_concordance_threshold
  }

  if (size(IlluminaGenotypingArray.gtc) == 0) {
    call InternalTasks.GenerateEmptyVariantCallingMetricsFile {
      input:
        chip_well_barcode = chip_well_barcode,
        sample_alias = sample_alias,
        chip_type = chip_type,
        reported_gender = reported_gender,
        autocall_version = autocall_version,
        output_metrics_basename = sample_alias,
        cluster_filename = cluster_filename,
        analysis_version_number = analysis_version_number,
        preemptible_tries = preemptible_tries
    }

    call InternalTasks.UploadArraysMetrics as UploadEmptyArraysMetrics {
      input:
        arrays_variant_calling_detail_metrics = GenerateEmptyVariantCallingMetricsFile.detail_metrics,
        disk_size = disk_size,
        preemptible_tries = preemptible_tries,
        authentication = authentication_block,
        service_account_filename = service_account_filename
    }

    call InternalTasks.BlacklistBarcode as BlacklistFailedNormalization {
      input:
        upload_metrics_output = UploadEmptyArraysMetrics.upload_metrics_empty_file,
        analysis_version = analysis_version_number,
        chip_well_barcode = chip_well_barcode,
        preemptible_tries = preemptible_tries,
        reason = "DATA_QUALITY",
        notes = "Normalization Failed",
        authentication = authentication_block,
        service_account_filename = service_account_filename
    }
  }

  if (size(IlluminaGenotypingArray.gtc) > 0) {
    call InternalTasks.VcfToMercuryFingerprintJson {
      input:
        input_vcf_file = select_first([IlluminaGenotypingArray.output_fingerprint_vcf]),
        input_vcf_index_file = select_first([IlluminaGenotypingArray.output_fingerprint_vcf_index]),
        variant_calling_detail_metrics_file = select_first([IlluminaGenotypingArray.arrays_variant_calling_detail_metrics]),
        sample_lsid = sample_lsid,
        output_json_filename = chip_well_barcode + ".fingerprint.json",
        disk_size = disk_size,
        preemptible_tries = preemptible_tries
    }

    if (defined(IlluminaGenotypingArray.bafregress_results_file)) {
      call InternalTasks.CreateBafRegressMetricsFile {
        input:
          input_file = select_first([IlluminaGenotypingArray.bafregress_results_file]),
          output_metrics_basefilename = chip_well_barcode,
          disk_size = disk_size,
          preemptible_tries = preemptible_tries
      }
    }

    call InternalTasks.UploadArraysMetrics {
      input:
        arrays_variant_calling_detail_metrics = select_first([IlluminaGenotypingArray.arrays_variant_calling_detail_metrics]),
        arrays_variant_calling_summary_metrics = IlluminaGenotypingArray.arrays_variant_calling_summary_metrics,
        arrays_control_code_summary_metrics = IlluminaGenotypingArray.arrays_variant_calling_control_metrics,
        fingerprinting_detail_metrics = IlluminaGenotypingArray.fingerprint_detail_metrics,
        fingerprinting_summary_metrics = IlluminaGenotypingArray.fingerprint_summary_metrics,
        genotype_concordance_summary_metrics = IlluminaGenotypingArray.genotype_concordance_summary_metrics,
        genotype_concordance_detail_metrics  = IlluminaGenotypingArray.genotype_concordance_detail_metrics,
        genotype_concordance_contingency_metrics = IlluminaGenotypingArray.genotype_concordance_contingency_metrics,
        verify_id_metrics = IlluminaGenotypingArray.contamination_metrics,
        bafregress_metrics = CreateBafRegressMetricsFile.output_metrics_file,
        disk_size = disk_size,
        preemptible_tries = preemptible_tries,
        authentication = authentication_block,
        service_account_filename = service_account_filename
    }

    if (select_first([IlluminaGenotypingArray.genotype_concordance_failed, false])) {
      call InternalTasks.BlacklistBarcode as BlacklistFailedGenotypeConcordance {
        input:
          upload_metrics_output = UploadArraysMetrics.upload_metrics_empty_file,
          analysis_version = analysis_version_number,
          chip_well_barcode = chip_well_barcode,
          preemptible_tries = preemptible_tries,
          reason = "GENOTYPE_CONCORDANCE",
          notes = "Genotype concordance below threshold: ~{genotype_concordance_threshold}",
          authentication = authentication_block,
          service_account_filename = service_account_filename
      }
    }

    if (defined(fingerprint_genotypes_vcf_file) &&
        select_first([IlluminaGenotypingArray.check_fingerprint_lod, 0.0]) < -3.0) {
      call InternalTasks.BlacklistBarcode as BlacklistFailedFingerprint {
        input:
          upload_metrics_output = UploadArraysMetrics.upload_metrics_empty_file,
          analysis_version = analysis_version_number,
          chip_well_barcode = chip_well_barcode,
          preemptible_tries = preemptible_tries,
          reason = "SAMPLE_MIXUP",
          notes = "Fingerprint LOD below -3.0",
          authentication = authentication_block,
          service_account_filename = service_account_filename
      }
    }
  }

  output {
    File GtcFile = IlluminaGenotypingArray.gtc
    File RedIdatMd5CloudPath = IlluminaGenotypingArray.red_idat_md5_cloud_path
    File GreenIdatMd5CloudPath = IlluminaGenotypingArray.green_idat_md5_cloud_path
    File? OutputVcfMd5CloudPath = IlluminaGenotypingArray.output_vcf_md5_cloud_path
    File? OutputVcfFile = IlluminaGenotypingArray.output_vcf
    File? OutputVcfIndexFile = IlluminaGenotypingArray.output_vcf_index
    File? BafRegressMetricsFile = CreateBafRegressMetricsFile.output_metrics_file
    File? ContaminationMetricsFile = IlluminaGenotypingArray.contamination_metrics
    File? ReferenceFingerprintVcf = DownloadGenotypes.output_vcf
    File? ReferenceFingerprintVcfIndex = DownloadGenotypes.output_vcf_index
    File? OutputFingerprintVcfFile = IlluminaGenotypingArray.output_fingerprint_vcf
    File? OutputFingerprintVcfIndexFile = IlluminaGenotypingArray.output_fingerprint_vcf_index
    File? OutputFingerprintJsonFile = VcfToMercuryFingerprintJson.output_json_file
    File ArraysVariantCallingDetailMetricsFile = select_first([IlluminaGenotypingArray.arrays_variant_calling_detail_metrics, GenerateEmptyVariantCallingMetricsFile.detail_metrics])
    File? ArraysVariantCallingSummaryMetricsFile = IlluminaGenotypingArray.arrays_variant_calling_summary_metrics
    File? ArraysVariantCallingControlMetricsFile = IlluminaGenotypingArray.arrays_variant_calling_control_metrics
    File? ArraysSubsetVariantCallingDetailMetricsFile = IlluminaGenotypingArray.arrays_subset_variant_calling_detail_metrics
    File? ArraysSubsetVariantCallingSummaryMetricsFile = IlluminaGenotypingArray.arrays_subset_variant_calling_summary_metrics
    File? ArraysSubsetVariantCallingControlMetricsFile = IlluminaGenotypingArray.arrays_subset_variant_calling_control_metrics
    File? FingerprintDetailMetricsFile = IlluminaGenotypingArray.fingerprint_detail_metrics
    File? FingerprintSummaryMetricsFile = IlluminaGenotypingArray.fingerprint_summary_metrics
    File? GenotypeConcordanceSummaryMetricsFile = IlluminaGenotypingArray.genotype_concordance_summary_metrics
    File? GenotypeConcordanceDetailMetricsFile  = IlluminaGenotypingArray.genotype_concordance_detail_metrics
    File? GenotypeConcordanceContingencyMetricsFile = IlluminaGenotypingArray.genotype_concordance_contingency_metrics
  }
  meta {
    allowNestedInputs: true
  }
}

task DownloadGenotypes {
  input {
    String sample_alias
    String sample_lsid
    String output_vcf_name

    File haplotype_database_file

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    Array[String] source_block

    Int preemptible_tries
  }

  String base_vcf = basename(output_vcf_name)
  Boolean is_compressed = basename(base_vcf, "gz") != base_vcf
  String vcf_suffix = if is_compressed then ".vcf.gz" else ".vcf"
  String vcf_index_suffix = if is_compressed then ".tbi" else ".idx"
  String output_vcf = base_vcf + vcf_suffix
  String output_vcf_index = output_vcf + vcf_index_suffix


  command <<<
    set -eo pipefail

    # TODO - using this won't work until we give the appropriate permissions to this user.
    AUTH=~{write_lines(source_block)} && source $AUTH

    # TODO - need to pull the credentials from the VAULT - but the permissions aren't working!
    #    --CREDENTIALS_VAULT_PATH secret/dsde/gotc/dev/wdl/secrets \

    java -Xms2g -Dpicard.useLegacyParser=false -jar /usr/gitc/picard-private.jar \
    DownloadGenotypes \
    --SAMPLE_ALIAS "~{sample_alias}" \
    --SAMPLE_LSID ~{sample_lsid} \
    --OUTPUT ~{output_vcf} \
    --CREATE_INDEX true \
    --REFERENCE_SEQUENCE ~{ref_fasta} \
    --HAPLOTYPE_MAP ~{haplotype_database_file} \
    # TODO - this should only be set to FLUIDIGM if it's 'clinical'
    # I think the tool may be broken if you are supplying nothing.  Need to test.
    --EXPECTED_GENOTYPING_PLATFORMS FLUIDIGM \
    --IGNORE_SPECIFIC_GENOTYPES_PLATFORM GENERAL_ARRAY \
    --IGNORE_SPECIFIC_GENOTYPES_LSID ~{sample_lsid} \
    --MERCURY_FP_STORE_URI https://portals.broadinstitute.org/portal-test/mercury-ws/fingerprint \
    --USERNAME dsp-fp-dev \
    --PASSWORD 'peph#datHe!2'

    # TODO - if it doesn't find any fingerprints, it generates no output (I think there's a different return code).
    # Need to then?? make an empty VCF and index and pass them along through subsequent steps - or have another flag
    # that is, no fingerprint, don't do VCFToFingerprintJson and don't do upload fingerprint.
    # so, really optional outputs I guess?
    # and you should add a test case for that.
  >>>

  runtime {
    docker: "us.gcr.io/broad-arrays-prod/arrays-picard-private:4.0.10-1602016912"
    memory: "3.5 GiB"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf = "~{output_vcf}"
    File output_vcf_index = "~{output_vcf_index}"
  }
}
