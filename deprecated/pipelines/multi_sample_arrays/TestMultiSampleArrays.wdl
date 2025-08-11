version 1.0


import "../../pipelines/wdl/arrays/multi_sample/MultiSampleArrays.wdl" as MultiSampleArrays
import "../../verification/VerifyMultiSampleArrays.wdl" as VerifyMultiSampleArrays
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/CopyFilesFromCloudToCloud.wdl" as Copy

workflow TestMultiSampleArrays {

    input {
      File samples_fofn
      File sample_indices_fofn
      File ref_fasta
      File ref_fasta_index
      File ref_dict
      String callset_name
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
  
    call MultiSampleArrays.MultiSampleArrays {
      input:
        samples_fofn = samples_fofn,
        sample_indices_fofn = sample_indices_fofn,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        callset_name = callset_name,
        preemptible_tries = preemptible_tries
  
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    MultiSampleArrays.combined_vcf_index,
                                    MultiSampleArrays.combined_vcf,
                                    ],
                                    
    ])

    

    # Copy results of pipeline to test results bucket
    call Copy.CopyFilesFromCloudToCloud as CopyToTestResults {
      input:
        files_to_copy             = flatten([pipeline_outputs]),
        vault_token_path          = vault_token_path,
        google_account_vault_path = google_account_vault_path,
        destination_cloud_path    = results_path
    }
  
    # If updating truth then copy output to truth bucket
    if (update_truth){
      call Copy.CopyFilesFromCloudToCloud as CopyToTruth {
        input: 
          files_to_copy             = flatten([pipeline_outputs]),
          vault_token_path          = vault_token_path,
          google_account_vault_path = google_account_vault_path,
          destination_cloud_path    = truth_path
      }
    }

    # If not updating truth then we need to collect all input for the validation WDL
    # This is achieved by passing each desired file/array[files] to GetValidationInputs
    if (!update_truth){
        call Utilities.GetValidationInputs as GetVcf {
          input:
            input_file = MultiSampleArrays.combined_vcf,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifyMultiSampleArrays.VerifyMultiSampleArrays as Verify {
        input:
          truth_vcf = GetVcf.truth_file, 
          test_vcf = GetVcf.results_file,
          done = CopyToTestResults.done
      }
    }
}