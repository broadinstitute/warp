version 1.0


import "../../pipelines/skylab/snM3C/snM3C.wdl" as snM3C
import "../../verification/VerifysnM3C.wdl" as VerifysnM3C
import "../../tasks/broad/Utilities.wdl" as Utilities
import "../../tasks/broad/CopyFilesFromCloudToCloud.wdl" as Copy

workflow TestsnM3C {

    input {
      Array[File] fastq_input_read1
      Array[File] fastq_input_read2
      File random_primer_indexes
      String plate_id
      String output_basename = plate_id
      File tarred_index_files
      File mapping_yaml
      File snakefile
      File chromosome_sizes
      File genome_fa

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
  
    call snM3C.snM3C {
      input:
        fastq_input_read1 = fastq_input_read1,
        fastq_input_read2 = fastq_input_read2,
        random_primer_indexes = random_primer_indexes,
        plate_id = plate_id,
        output_basename = output_basename,
        tarred_index_files = tarred_index_files,
        mapping_yaml = mapping_yaml,
        snakefile = snakefile,
        chromosome_sizes = chromosome_sizes,
        genome_fa = genome_fa
  
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    snM3C.hicFiles,
                                    snM3C.detail_statsFiles,
                                    snM3C.bamFiles,
                                    snM3C.allc_CGNFiles,
                                    snM3C.allcFiles,
                                    snM3C.MappingSummary,
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

    # This is achieved by passing each desired file/array[files] to GetValidationInputs
    if (!update_truth){
        call Utilities.GetValidationInputs as GetMappingSummary {
          input:
            input_file = snM3C.MappingSummary,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifysnM3C.VerifysnM3C as Verify {
        input:
          truth_mapping_summary = GetMappingSummary.truth_file, 
          test_mapping_summary = GetMappingSummary.results_file,
          done = CopyToTestResults.done
      }
    }
}