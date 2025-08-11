version 1.0


import "../../pipelines/wdl/reprocessing/cram_to_unmapped_bams/CramToUnmappedBams.wdl" as CramToUnmappedBams
import "../../verification/VerifyCramToUnmappedBamsUpdated.wdl" as VerifyCramToUnmappedBamsUpdated
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestCramToUnmappedBams {

    input {
      File? input_cram
      File? input_bam
      File? ref_fasta
      File? ref_fasta_index
      File? output_map
      String base_file_name
      String unmapped_bam_suffix = ".unmapped.bam"
      Int additional_disk = 20

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
    }

    meta {
      allowNestedInputs: true
    }
  
    call CramToUnmappedBams.CramToUnmappedBams {
      input:
        input_cram = input_cram,
        input_bam = input_bam,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        output_map = output_map,
        base_file_name = base_file_name,
        unmapped_bam_suffix = unmapped_bam_suffix,
        additional_disk = additional_disk
  
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    # Array[File] outputs
                                    CramToUnmappedBams.unmapped_bams,
                                    CramToUnmappedBams.validation_report,
                                    
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

    # If not updating truth then we need to collect all input for the validation WDL
    # This is achieved by passing each desired file/array[files] to GetValidationInputs
    if (!update_truth){
        call Utilities.GetValidationInputs as GetBam {
          input:
            input_files = CramToUnmappedBams.unmapped_bams,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifyCramToUnmappedBamsUpdated.VerifyCramToUnmappedBams as Verify {
        input:
          truth_bam = GetBam.truth_files, 
          test_bam = GetBam.results_files,
          done = CopyToTestResults.done
      }
  
    }





}