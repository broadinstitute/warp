version 1.0


import "../../pipelines/wdl/qc/CheckFingerprint.wdl" as CheckFingerprint
import "../../verification/VerifyCheckFingerprint.wdl" as VerifyCheckFingerprint
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/CopyFilesFromCloudToCloud.wdl" as Copy

workflow TestCheckFingerprint {

    input {
      File? input_vcf
      File? input_vcf_index
      File? input_bam
      File? input_bam_index
      String? input_sample_alias
      Boolean read_fingerprint_from_mercury = false
      File? fingerprint_genotypes_vcf
      File? fingerprint_genotypes_vcf_index
      String? sample_lsid
      String sample_alias
      String output_basename
      File ref_fasta
      File ref_fasta_index
      File ref_dict
      File haplotype_database_file

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
      String vault_token_path
      String? vault_token_path_arrays
      String environment
      String google_account_vault_path
    }

    meta {
      allowNestedInputs: true
    }
  
    call CheckFingerprint.CheckFingerprint {
      input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        input_bam = input_bam,
        input_bam_index = input_bam_index,
        input_sample_alias = input_sample_alias,
        read_fingerprint_from_mercury = read_fingerprint_from_mercury,
        fingerprint_genotypes_vcf = fingerprint_genotypes_vcf,
        fingerprint_genotypes_vcf_index = fingerprint_genotypes_vcf_index,
        sample_lsid = sample_lsid,
        sample_alias = sample_alias,
        output_basename = output_basename,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        haplotype_database_file = haplotype_database_file,
        environment = environment,
        vault_token_path = vault_token_path
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    # File? outputs
                                    select_all([CheckFingerprint.reference_fingerprint_vcf_index]),
                                    select_all([CheckFingerprint.reference_fingerprint_vcf]),
                                    
    ])

    
    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    # File? outputs
                                    select_all([CheckFingerprint.fingerprint_detail_metrics_file]),
                                    select_all([CheckFingerprint.fingerprint_summary_metrics_file]),
                                    
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
        call Utilities.GetValidationInputs as GetMetrics {
          input:
            input_files = pipeline_metrics,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetFingerprintVcf {
          input:
            input_file = CheckFingerprint.reference_fingerprint_vcf,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifyCheckFingerprint.VerifyCheckFingerprint as Verify {
        input:
          truth_metrics = GetMetrics.truth_files, 
          test_metrics = GetMetrics.results_files,
          truth_fingerprint_vcf = GetFingerprintVcf.truth_file, 
          test_fingerprint_vcf = GetFingerprintVcf.results_file,
          done = CopyToTestResults.done
      }
    }
}