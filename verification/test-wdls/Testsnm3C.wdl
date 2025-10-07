version 1.0


import "../../pipelines/wdl/snm3C/snm3C.wdl" as snm3C
import "../../verification/Verifysnm3C.wdl" as Verifysnm3C
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow Testsnm3C {

    input {
      Array[File] fastq_input_read1
      Array[File] fastq_input_read2
      File random_primer_indexes
      String plate_id
      String cloud_provider
      File tarred_index_files
      File genome_fa
      File chromosome_sizes
      String r1_adapter = "AGATCGGAAGAGCACACGTCTGAAC"
      String r2_adapter = "AGATCGGAAGAGCGTCGTGTAGGGA"
      Int r1_left_cut = 10
      Int r1_right_cut = 10
      Int r2_left_cut = 10
      Int r2_right_cut = 10
      Int min_read_length = 30
      Int num_upstr_bases = 0
      Int num_downstr_bases = 2
      Int compress_level = 5
      Int batch_number = 6

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
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
          cloud_provider = cloud_provider,
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
          batch_number = batch_number

    }


    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    snm3C.MappingSummary,
                                    ],
                                    # Array[File] outputs
                                    snm3C.reference_version,
                                    snm3C.unique_reads_cgn_extraction_allc_array,
                                    snm3C.unique_reads_cgn_extraction_tbi_array,
                                    snm3C.unique_reads_cgn_extraction_allc_extract_array,
                                    snm3C.unique_reads_cgn_extraction_tbi_extract_array,
                                    snm3C.name_sorted_bam_array,
                                    snm3C.all_reads_3C_contacts_array
    ])

    

    # Copy results of pipeline to test results bucket
    call Copy.TerraCopyFilesFromCloudToCloud as CopyToTestResults {
      input:
        files_to_copy             = flatten([pipeline_outputs]),
        destination_cloud_path    = results_path
    }
  
    # If updating truth then copy output to truth bucket
    if (update_truth){
      call Copy.TerraCopyFilesFromCloudToCloud as CopyToTruth {
        input: 
          files_to_copy             = flatten([pipeline_outputs]),
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
        call Utilities.GetValidationInputs as GetNameSortedBamArray {
            input:
                input_files = snm3C.name_sorted_bam_array,
                results_path = results_path,
                truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetUniqueReadsCgnExtractionAllcArray {
            input:
                input_files = snm3C.unique_reads_cgn_extraction_allc_array,
                results_path = results_path,
                truth_path = truth_path
        }

        call Utilities.GetValidationInputs as GetUniqueReadsCgnExtractionAllcExtractArray {
            input:
                input_files = snm3C.unique_reads_cgn_extraction_allc_extract_array,
                results_path = results_path,
                truth_path = truth_path
        }

        call Utilities.GetValidationInputs as GetAllReads3CContactsArray {
            input:
                input_files = snm3C.all_reads_3C_contacts_array,
                results_path = results_path,
                truth_path = truth_path
        }

        call Verifysnm3C.Verifysnm3C as Verify {
            input:
                truth_mapping_summary = GetMappingSummary.truth_file,
                test_mapping_summary = GetMappingSummary.results_file,
                truth_name_sorted_bam_array = GetNameSortedBamArray.truth_files,
                test_name_sorted_bam_array = GetNameSortedBamArray.results_files,
                truth_unique_reads_cgn_extraction_allc_array = GetUniqueReadsCgnExtractionAllcArray.truth_files,
                test_unique_reads_cgn_extraction_allc_array = GetUniqueReadsCgnExtractionAllcArray.results_files,
                truth_unique_reads_cgn_extraction_allc_extract_array = GetUniqueReadsCgnExtractionAllcExtractArray.truth_files,
                test_unique_reads_cgn_extraction_allc_extract_array = GetUniqueReadsCgnExtractionAllcExtractArray.results_files,
                truth_all_reads_3C_contacts_array = GetAllReads3CContactsArray.truth_files,
                test_all_reads_3C_contacts_array = GetAllReads3CContactsArray.results_files,
                done = CopyToTestResults.done
        }
    }
}