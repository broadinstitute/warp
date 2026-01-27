version 1.0

import "../../../../tasks/wdl/IlluminaGenotypingArrayTasks.wdl" as GenotypingTasks
import "../../../../tasks/wdl/Qc.wdl" as Qc

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

workflow IlluminaGenotypingArray {

  String pipeline_version = "1.12.27"

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

    File extended_chip_manifest_file
    File cluster_file
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

    Float genotype_concordance_threshold = 0.95
  }

  call GenotypingTasks.AutoCall {
    input:
      chip_well_barcode = chip_well_barcode,
      green_idat_cloud_path = green_idat_cloud_path,
      red_idat_cloud_path = red_idat_cloud_path,
      bead_pool_manifest_file = bead_pool_manifest_file,
      cluster_file = cluster_file,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
  }

  call GenotypingTasks.Md5Sum as RedIdatMd5Sum {
    input:
      input_cloud_path = red_idat_cloud_path
  }

  call GenotypingTasks.Md5Sum as GreenIdatMd5Sum {
    input:
      input_cloud_path = green_idat_cloud_path
  }

  if (defined(gender_cluster_file)) {
    call GenotypingTasks.AutoCall as GenderAutocall {
      input:
        chip_well_barcode = chip_well_barcode,
        green_idat_cloud_path = green_idat_cloud_path,
        red_idat_cloud_path = red_idat_cloud_path,
        bead_pool_manifest_file = bead_pool_manifest_file,
        cluster_file = gender_cluster_file,
        is_gender_autocall = true,
        disk_size = disk_size,
        preemptible_tries = preemptible_tries
    }
  }

  if (size(AutoCall.gtc_file) > 0) {

    call GenotypingTasks.GtcToVcf {
      input:
        vcf_filename = chip_well_barcode + ".vcf.gz",
        input_gtc = AutoCall.gtc_file,
        gender_gtc = GenderAutocall.gtc_file,
        extended_chip_manifest_file = extended_chip_manifest_file,
        cluster_file = cluster_file,
        bead_pool_manifest_file = bead_pool_manifest_file,
        sample_alias = sample_alias,
        analysis_version_number = analysis_version_number,
        reported_gender = reported_gender,
        fingerprint_genotypes_vcf_file = fingerprint_genotypes_vcf_file,
        fingerprint_genotypes_vcf_index_file = fingerprint_genotypes_vcf_index_file,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        disk_size = disk_size,
        preemptible_tries = preemptible_tries,
        pipeline_version = "IlluminaGenotypingArray_v" + pipeline_version
    }

    if (defined(minor_allele_frequency_file)) {
      call GenotypingTasks.BafRegress {
        input:
          input_vcf = GtcToVcf.output_vcf,
          input_vcf_index = GtcToVcf.output_vcf_index,
          maf_file = minor_allele_frequency_file,
          output_results_filename = chip_well_barcode + ".results.txt",
          disk_size = disk_size,
          preemptible_tries = preemptible_tries,
      }
    }

    call GenotypingTasks.VcfToAdpc {
      input:
        input_vcf = GtcToVcf.output_vcf,
        contamination_controls_vcf = contamination_controls_vcf,
        output_adpc_filename = chip_well_barcode + ".adpc.bin",
        samples_filename = chip_well_barcode + ".samples.txt",
        num_markers_filename = chip_well_barcode + ".num_markers.txt",
        disk_size = disk_size,
        preemptible_tries = preemptible_tries
    }

    call GenotypingTasks.VerifyIDIntensity {
      input:
        input_vcf = GtcToVcf.output_vcf,
        input_adpc_file = VcfToAdpc.output_adpc,
        num_samples = length(read_lines(VcfToAdpc.samples_file)),
        num_markers = VcfToAdpc.num_markers,
        output_filename = chip_well_barcode + ".verifyIDIntensity.txt",
        disk_size = disk_size,
        preemptible_tries = preemptible_tries
    }

    call GenotypingTasks.CreateVerifyIDIntensityContaminationMetricsFile {
      input:
        expected_first_sample_name = chip_well_barcode,
        samples_file = VcfToAdpc.samples_file,
        input_file = VerifyIDIntensity.output_file,
        output_metrics_basefilename = chip_well_barcode,
        disk_size = disk_size,
        preemptible_tries = preemptible_tries
    }

    if (defined(zcall_thresholds_file)) {
      call GenotypingTasks.BpmToNormalizationManifestCsv {
        input:
          bead_pool_manifest_file = bead_pool_manifest_file,
          cluster_file = cluster_file,
          bead_pool_manifest_csv_file = chip_type + ".bpm.csv",
          disk_size = disk_size,
          preemptible_tries = preemptible_tries
      }

      call GenotypingTasks.zCall {
        input:
          zcall_ped_filename = chip_well_barcode + ".zcall.ped",
          zcall_map_filename = chip_well_barcode + ".zcall.map",
          input_gtc = AutoCall.gtc_file,
          bead_pool_manifest_csv_file = BpmToNormalizationManifestCsv.output_file,
          zcall_thresholds_file = select_first([zcall_thresholds_file]),
          disk_size = disk_size,
          preemptible_tries = preemptible_tries
      }

      call GenotypingTasks.MergePedIntoVcf {
        input:
          input_vcf = GtcToVcf.output_vcf,
          input_vcf_index = GtcToVcf.output_vcf_index,
          output_vcf_filename = sub(GtcToVcf.output_vcf, "gs://.*/", ""),
          ped_file = zCall.ped_file,
          map_file = zCall.map_file,
          zcall_thresholds_file = select_first([zcall_thresholds_file]),
          zcall_version = "1.0.0.0",
          disk_size = disk_size,
          preemptible_tries = preemptible_tries
      }
    }

    # if zCall doesn't run, then MergePedIntoVcf doesn't run.
    # if MergePedIntoVcf doesn't run, then MergePedIntoVcf.output_vcf and output_vcf_index should just be
    # the MergePedIntoVcf.input_vcf and input_vcf_index
    # using select_first to cast File? to File type
    File final_output_vcf =       select_first([MergePedIntoVcf.output_vcf, GtcToVcf.output_vcf])
    File final_output_vcf_index = select_first([MergePedIntoVcf.output_vcf_index, GtcToVcf.output_vcf_index])

    call GenotypingTasks.Md5Sum as VcfMd5Sum {
      input:
        input_cloud_path = final_output_vcf
    }

    call GenotypingTasks.CollectArraysVariantCallingMetrics {
      input:
        input_vcf_file = final_output_vcf,
        input_vcf_index_file = final_output_vcf_index,
        dbSNP_vcf_file = dbSNP_vcf,
        dbSNP_vcf_index_file = dbSNP_vcf_index,
        call_rate_threshold = call_rate_threshold,
        output_metrics_basename = chip_well_barcode,
        disk_size = disk_size,
        preemptible_tries = preemptible_tries
    }

    if (defined(subsampled_metrics_interval_list)) {
      call GenotypingTasks.SubsetArrayVCF {
        input:
          intervals = select_first([subsampled_metrics_interval_list]),
          input_vcf_file = final_output_vcf,
          input_vcf_index_file = final_output_vcf_index,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict
      }

      call GenotypingTasks.CollectArraysVariantCallingMetrics as SubsettedVariantCallingMetrics {
        input:
          input_vcf_file = SubsetArrayVCF.output_vcf,
          input_vcf_index_file = SubsetArrayVCF.output_vcf_index,
          dbSNP_vcf_file = dbSNP_vcf,
          dbSNP_vcf_index_file = dbSNP_vcf_index,
          call_rate_threshold = call_rate_threshold,
          output_metrics_basename = chip_well_barcode + "_subset",
          disk_size = disk_size,
          preemptible_tries = preemptible_tries
        }
    }
    call GenotypingTasks.SelectVariants as SelectFingerprintVariants {
      input:
        input_vcf_file = final_output_vcf,
        input_vcf_index_file = final_output_vcf_index,
        variant_rsids_file = variant_rsids_file,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        output_vcf_filename = chip_well_barcode + ".fingerprint.vcf.gz",
        preemptible_tries = preemptible_tries
    }
    if (defined(fingerprint_genotypes_vcf_file)) {
      call Qc.CheckFingerprintTask as CheckFingerprintTask{
        input:
          input_vcf = final_output_vcf,
          input_vcf_index = final_output_vcf_index,
          input_sample_alias = chip_well_barcode,
          genotypes = select_first([fingerprint_genotypes_vcf_file]),
          genotypes_index = select_first([fingerprint_genotypes_vcf_index_file]),
          expected_sample_alias = sample_alias,
          output_basename = chip_well_barcode,
          genotype_lod_threshold = 1.9,
              # Paraphrased from Yossi:
              # Override the default LOD threshold of 5 because if the PL field
              # is missing from the VCF, CheckFingerprint will default to an error
              # rate equivalent to a LOD score of 2, and we don't want to see
              # confident LOD scores w/ no confident SNPs.
          haplotype_database_file = haplotype_database_file,
          preemptible_tries = preemptible_tries
      }
    }

    if (defined(control_sample_vcf_file) && defined(control_sample_intervals_file) && defined(control_sample_name)) {
      call GenotypingTasks.VcfToIntervalList {
        input:
          vcf_file = GtcToVcf.output_vcf,
          interval_list_filename = chip_well_barcode + ".interval_list",
          disk_size = disk_size,
          preemptible_tries = preemptible_tries
      }

      call GenotypingTasks.SelectVariants as SelectVariantsForGenotypeConcordance {
        input:
          input_vcf_file = GtcToVcf.output_vcf,
          input_vcf_index_file = GtcToVcf.output_vcf_index,
          output_vcf_filename = chip_well_barcode + ".select_variants.vcf",
          excludeFiltered = true,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          preemptible_tries = preemptible_tries
      }

      call GenotypingTasks.GenotypeConcordance {
        input:
          call_vcf_file = SelectVariantsForGenotypeConcordance.output_vcf,
          call_vcf_index_file = SelectVariantsForGenotypeConcordance.output_vcf_index,
          call_intervals_file = VcfToIntervalList.interval_list_file,
          call_sample_name = chip_well_barcode,
          truth_vcf_file = control_sample_vcf_file,
          truth_vcf_index_file = control_sample_vcf_index_file,
          truth_intervals_file = control_sample_intervals_file,
          truth_sample_name = control_sample_name,
          genotype_concordance_threshold = genotype_concordance_threshold,
          output_metrics_basename = chip_well_barcode,
          disk_size = disk_size,
          preemptible_tries = preemptible_tries
      }
    }
  }

  output {
    String chip_well_barcode_output = chip_well_barcode
    Int analysis_version_number_output = analysis_version_number
    String autocall_version = AutoCall.autocall_version
    File gtc = AutoCall.gtc_file
    File red_idat_md5_cloud_path = RedIdatMd5Sum.md5_cloud_path
    File green_idat_md5_cloud_path = GreenIdatMd5Sum.md5_cloud_path
    File? output_vcf_md5_cloud_path = VcfMd5Sum.md5_cloud_path
    File? output_vcf = final_output_vcf
    File? output_vcf_index = final_output_vcf_index
    File? bafregress_results_file = BafRegress.results_file
    File? contamination_metrics = CreateVerifyIDIntensityContaminationMetricsFile.output_metrics_file
    File? output_fingerprint_vcf = SelectFingerprintVariants.output_vcf
    File? output_fingerprint_vcf_index = SelectFingerprintVariants.output_vcf_index
    File? arrays_variant_calling_detail_metrics = CollectArraysVariantCallingMetrics.detail_metrics
    File? arrays_variant_calling_summary_metrics = CollectArraysVariantCallingMetrics.summary_metrics
    File? arrays_variant_calling_control_metrics = CollectArraysVariantCallingMetrics.control_metrics
    File? arrays_subset_variant_calling_detail_metrics = SubsettedVariantCallingMetrics.detail_metrics
    File? arrays_subset_variant_calling_summary_metrics = SubsettedVariantCallingMetrics.summary_metrics
    File? arrays_subset_variant_calling_control_metrics = SubsettedVariantCallingMetrics.control_metrics
    File? fingerprint_detail_metrics = CheckFingerprintTask.detail_metrics
    File? fingerprint_summary_metrics = CheckFingerprintTask.summary_metrics
    Float? check_fingerprint_lod = CheckFingerprintTask.lod
    File? genotype_concordance_summary_metrics = GenotypeConcordance.summary_metrics
    File? genotype_concordance_detail_metrics = GenotypeConcordance.detail_metrics
    File? genotype_concordance_contingency_metrics = GenotypeConcordance.contingency_metrics
    Boolean? genotype_concordance_failed = GenotypeConcordance.fails_concordance
  }
  meta {
    allowNestedInputs: true
  }
}
