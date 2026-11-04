version 1.0


import "../../pipelines/wdl/genotyping/illumina/IlluminaGenotypingArray.wdl" as IlluminaGenotypingArray
import "../../verification/VerifyIlluminaGenotypingArray.wdl" as VerifyIlluminaGenotypingArray
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestIlluminaGenotypingArray {

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
      File? fingerprint_genotypes_vcf_file
      File? fingerprint_genotypes_vcf_index_file
      File haplotype_database_file
      File variant_rsids_file
      File? subsampled_metrics_interval_list
      File? contamination_controls_vcf
      File? minor_allele_frequency_file
      File? control_sample_vcf_file
      File? control_sample_vcf_index_file
      File? control_sample_intervals_file
      String? control_sample_name
      Int disk_size
      Int preemptible_tries
      Float genotype_concordance_threshold = 0.95

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
    }

    meta {
      allowNestedInputs: true
    }


    call IlluminaGenotypingArray.IlluminaGenotypingArray {
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
        extended_chip_manifest_file = extended_chip_manifest_file,
        cluster_file = cluster_file,
        gender_cluster_file = gender_cluster_file,
        zcall_thresholds_file = zcall_thresholds_file,
        fingerprint_genotypes_vcf_file = fingerprint_genotypes_vcf_file,
        fingerprint_genotypes_vcf_index_file = fingerprint_genotypes_vcf_index_file,
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


    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    IlluminaGenotypingArray.green_idat_md5_cloud_path,
                                    IlluminaGenotypingArray.red_idat_md5_cloud_path,
                                    IlluminaGenotypingArray.gtc,
                                    ],
                                    # File? outputs
                                    select_all([IlluminaGenotypingArray.output_fingerprint_vcf_index]),
                                    select_all([IlluminaGenotypingArray.output_fingerprint_vcf]),
                                    select_all([IlluminaGenotypingArray.bafregress_results_file]),
                                    select_all([IlluminaGenotypingArray.output_vcf_index]),
                                    select_all([IlluminaGenotypingArray.output_vcf]),
                                    select_all([IlluminaGenotypingArray.output_vcf_md5_cloud_path]),
    ])


    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    # File? outputs
                                    select_all([IlluminaGenotypingArray.genotype_concordance_contingency_metrics]),
                                    select_all([IlluminaGenotypingArray.genotype_concordance_detail_metrics]),
                                    select_all([IlluminaGenotypingArray.genotype_concordance_summary_metrics]),
                                    select_all([IlluminaGenotypingArray.fingerprint_summary_metrics]),
                                    select_all([IlluminaGenotypingArray.fingerprint_detail_metrics]),
                                    select_all([IlluminaGenotypingArray.arrays_subset_variant_calling_control_metrics]),
                                    select_all([IlluminaGenotypingArray.arrays_subset_variant_calling_summary_metrics]),
                                    select_all([IlluminaGenotypingArray.arrays_subset_variant_calling_detail_metrics]),
                                    select_all([IlluminaGenotypingArray.arrays_variant_calling_control_metrics]),
                                    select_all([IlluminaGenotypingArray.arrays_variant_calling_summary_metrics]),
                                    select_all([IlluminaGenotypingArray.arrays_variant_calling_detail_metrics]),
                                    select_all([IlluminaGenotypingArray.contamination_metrics]),
    ])

    # Copy results of pipeline to test results bucket
    call Copy.TerraCopyFilesFromCloudToCloud as CopyToTestResults {
      input:
        files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
        destination_cloud_path    = results_path
    }

    # If updating truth then copy output to truth bucket
    if (update_truth){
      call Copy.TerraCopyFilesFromCloudToCloud as CopyToTruth {
        input:
          files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
          destination_cloud_path    = truth_path
      }
    }

    # If not updating truth then we need to collect all input for the validation WDL
    # This is achieved by passing each desired file/array[files] to GetValidationInputs
    if (!update_truth){
        call Utilities.GetValidationInputs as GetMetrics {
          input:
            input_files = pipeline_metrics,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetGtc {
          input:
            input_file = IlluminaGenotypingArray.gtc,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetVcf {
          input:
            input_file = select_first([IlluminaGenotypingArray.output_vcf,""]),
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetFpVcf {
          input:
            input_file = select_first([IlluminaGenotypingArray.output_fingerprint_vcf,""]),
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetRedIdatMd5 {
          input:
            input_file = IlluminaGenotypingArray.red_idat_md5_cloud_path,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetGreenIdatMd5 {
          input:
            input_file = IlluminaGenotypingArray.green_idat_md5_cloud_path,
            results_path = results_path,
            truth_path = truth_path
        }
      call VerifyIlluminaGenotypingArray.VerifyIlluminaGenotypingArray as Verify {
        input:
          truth_metrics = GetMetrics.truth_files,
          test_metrics = GetMetrics.results_files,
          truth_gtc = GetGtc.truth_file,
          test_gtc = GetGtc.results_file,
          truth_vcf = GetVcf.truth_file,
          test_vcf = GetVcf.results_file,
          truth_fp_vcf = GetFpVcf.truth_file,
          test_fp_vcf = GetFpVcf.results_file,
          truth_red_idat_md5 = GetRedIdatMd5.truth_file,
          test_red_idat_md5 = GetRedIdatMd5.results_file,
          truth_green_idat_md5 = GetGreenIdatMd5.truth_file,
          test_green_idat_md5 = GetGreenIdatMd5.results_file,
          bead_pool_manifest_file = bead_pool_manifest_file,
          done = CopyToTestResults.done
      }
    }

    output {
    }
}