version 1.0


import "../../pipelines/skylab/snm3C/snm3C.wdl" as snm3C
import "../../verification/Verifysnm3C.wdl" as Verifysnm3C
import "../../tasks/broad/Utilities.wdl" as Utilities
import "../../tasks/broad/CopyFilesFromCloudToCloud.wdl" as Copy

workflow Testsnm3C {

    input {
      Array[File] fastq_input_read1
      Array[File] fastq_input_read2
      File random_primer_indexes
      String plate_id
      File tarred_index_files
      File genome_fa
      File chromosome_sizes
      String r1_adapter = "AGATCGGAAGAGCACACGTCTGAAC"
      String r2_adapter = "AGATCGGAAGAGCGTCGTGTAGGGA"
      #Int batch_number
      Int r1_left_cut = 10
      Int r1_right_cut = 10
      Int r2_left_cut = 10
      Int r2_right_cut = 10
      Int min_read_length = 30
      Int num_upstr_bases = 0
      Int num_downstr_bases = 2
      Int compress_level = 5
      Int batch_number

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
      String vault_token_path
      String google_account_vault_path

      String docker = "us.gcr.io/broad-gotc-prod/m3c-yap-hisat:2.4"
    }

    meta {
      allowNestedInputs: true
    }
  
    call snm3C.snm3C {
      input:
        fastq_input_read1 = fastq_input_read1,
        fastq_input_read2 = fastq_input_read2,
        random_primer_indexes = random_primer_indexes,
        plate_id = plate_id,
        tarred_index_files = tarred_index_files,
        genome_fa = genome_fa,
        chromosome_sizes = chromosome_sizes,
        r1_adapter = r1_adapter,
        r2_adapter = r2_adapter,
        r1_left_cut = r1_left_cut,
        r1_right_cut = r1_right_cut,
        r2_left_cut = r2_left_cut,
        r2_right_cut = r2_right_cut,
        min_read_length = min_read_length,
        num_upstr_bases = num_upstr_bases,
        num_downstr_bases = num_downstr_bases,
        compress_level = compress_level,
        batch_number = batch_number,
        docker = docker
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    snm3C.MappingSummary,
                                    ],
                                    # Array[File] outputs
                                    snm3C.reference_version,
                                    snm3C.unique_reads_cgn_extraction_allc,
                                    snm3C.unique_reads_cgn_extraction_tbi,
                                    snm3C.unique_reads_cgn_extraction_allc_extract,
                                    snm3C.unique_reads_cgn_extraction_tbi_extract,
                                    snm3C.name_sorted_bams,
                                    snm3C.all_reads_dedup_contacts,
                                    snm3C.all_reads_3C_contacts,
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
            input_file = snm3C.MappingSummary,
            results_path = results_path,
            truth_path = truth_path
        }

      call Verifysnm3C.Verifysnm3C as Verify {
        input:
          truth_mapping_summary = GetMappingSummary.truth_file, 
          test_mapping_summary = GetMappingSummary.results_file,
          done = CopyToTestResults.done
      }
    }
}