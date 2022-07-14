version 1.0


import "../../pipelines/broad/dna_seq/germline/joint_genotyping/UltimaGenomics/UltimaGenomicsJointGenotyping.wdl" as UltimaGenomicsJointGenotyping
import "../../verification/VerifyUltimaGenomicsJointGenotyping.wdl" as VerifyUltimaGenomicsJointGenotyping
import "../../tasks/broad/Utilities.wdl" as Utilities
import "../../tasks/broad/CopyFilesFromCloudToCloud.wdl" as Copy

workflow TestUltimaGenomicsJointGenotyping {

    input {
      File unpadded_intervals_file
      String callset_name
      File sample_name_map
      File ref_fasta
      File ref_fasta_index
      File ref_dict
      File dbsnp_vcf
      File dbsnp_vcf_index
      Int small_disk = 100
      Int medium_disk = 200
      Int large_disk = 1000
      Int huge_disk = 2000
      File haplotype_database
      File eval_interval_list
      Float excess_het_threshold = 54.69
      Int? top_level_scatter_count
      Boolean? gather_vcfs
      Float unbounded_scatter_count_scale_factor = 0.15
      Boolean cross_check_fingerprints = true
      Boolean scatter_cross_check_fingerprints = false

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
  
    call UltimaGenomicsJointGenotyping.UltimaGenomicsJointGenotyping {
      input:
        unpadded_intervals_file = unpadded_intervals_file,
        callset_name = callset_name,
        sample_name_map = sample_name_map,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        small_disk = small_disk,
        medium_disk = medium_disk,
        large_disk = large_disk,
        huge_disk = huge_disk,
        haplotype_database = haplotype_database,
        eval_interval_list = eval_interval_list,
        excess_het_threshold = excess_het_threshold,
        top_level_scatter_count = top_level_scatter_count,
        gather_vcfs = gather_vcfs,
        unbounded_scatter_count_scale_factor = unbounded_scatter_count_scale_factor,
        cross_check_fingerprints = cross_check_fingerprints,
        scatter_cross_check_fingerprints = scatter_cross_check_fingerprints
  
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    UltimaGenomicsJointGenotyping.unfiltered_sites_only_vcf_index,
                                    UltimaGenomicsJointGenotyping.unfiltered_sites_only_vcf,
                                    UltimaGenomicsJointGenotyping.crosscheck_fingerprint_check,
                                    ],
                                    # Array[File] outputs
                                    UltimaGenomicsJointGenotyping.output_intervals,
                                    UltimaGenomicsJointGenotyping.unfiltered_output_vcf_indices,
                                    UltimaGenomicsJointGenotyping.unfiltered_output_vcfs,
                                    
    ])

    
    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    [ # File outputs
                                    UltimaGenomicsJointGenotyping.unfiltered_summary_metrics_file,
                                    UltimaGenomicsJointGenotyping.unfiltered_detail_metrics_file,
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

    # This is achieved by passing each desired file/array[files] to GetValidationInputs
    if (!update_truth){
        call Utilities.GetValidationInputs as GetVcfs {
          input:
            input_files = UltimaGenomicsJointGenotyping.unfiltered_output_vcfs,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetVcfIndexes {
          input:
            input_files = UltimaGenomicsJointGenotyping.unfiltered_output_vcf_indices,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetIntervals {
          input:
            input_files = UltimaGenomicsJointGenotyping.output_intervals,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetMetrics {
          input:
            input_files = pipeline_metrics,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetFingerprint {
          input:
            input_file = UltimaGenomicsJointGenotyping.crosscheck_fingerprint_check,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifyUltimaGenomicsJointGenotyping.VerifyJointGenotyping as Verify {
        input:
          truth_vcfs = GetVcfs.truth_files, 
          test_vcfs = GetVcfs.results_files,
          truth_vcf_indexes = GetVcfIndexes.truth_files, 
          test_vcf_indexes = GetVcfIndexes.results_files,
          truth_intervals = GetIntervals.truth_files, 
          test_intervals = GetIntervals.results_files,
          truth_metrics = GetMetrics.truth_files, 
          test_metrics = GetMetrics.results_files,
          truth_fingerprint = GetFingerprint.truth_file, 
          test_fingerprint = GetFingerprint.results_file,
          done = CopyToTestResults.done
      }
    }
}