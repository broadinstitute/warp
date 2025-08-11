version 1.0


import "../../pipelines/wdl/scATAC/scATAC.wdl" as scATAC
import "../../verification/VerifyscATAC.wdl" as VerifyscATAC
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/CopyFilesFromCloudToCloud.wdl" as Copy

workflow TestscATAC {

    input {
      File input_fastq1
      File input_fastq2
      String input_id
      String genome_name
      File input_reference
      String output_bam = input_id + "_aligned.bam"
      String bin_size_list = "10000"

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
  
    call scATAC.scATAC {
      input:
        input_fastq1 = input_fastq1,
        input_fastq2 = input_fastq2,
        input_id = input_id,
        genome_name = genome_name,
        input_reference = input_reference,
        output_bam = output_bam,
        bin_size_list = bin_size_list
  
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    scATAC.breakout_barcodesSection,
                                    scATAC.breakout_binCounts,
                                    scATAC.breakout_binCoordinates,
                                    scATAC.breakout_fragments,
                                    scATAC.breakout_barcodes,
                                    scATAC.output_aligned_bam,
                                    scATAC.output_snap,
                                    scATAC.output_snap_qc,
                                    ],
                                    
    ])

    # Collect the matrix csv files
    Array[String] matrix_files = [
                                    scATAC.breakout_barcodesSection,
                                    scATAC.breakout_binCounts,
                                    scATAC.breakout_binCoordinates,
                                    scATAC.breakout_fragments,
                                    scATAC.breakout_barcodes,
                                  ]

    

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

    # This is achieved by passing each desired file/array[files] to GetValidationInputs
    if (!update_truth){
        call Utilities.GetValidationInputs as GetBam {
          input:
            input_file = scATAC.output_aligned_bam,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetMatrixFiles {
          input:
            input_files = matrix_files,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifyscATAC.VerifyscATAC as Verify {
        input:
          truth_bam = GetBam.truth_file, 
          test_bam = GetBam.results_file,
          truth_matrix_files = GetMatrixFiles.truth_files, 
          test_matrix_files = GetMatrixFiles.results_files,
          done = CopyToTestResults.done
      }
    }
}