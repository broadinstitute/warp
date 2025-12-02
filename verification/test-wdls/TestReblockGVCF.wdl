version 1.0


import "../../pipelines/wdl/dna_seq/germline/joint_genotyping/reblocking/ReblockGVCF.wdl" as ReblockGVCF
import "../../verification/VerifyGvcf.wdl" as VerifyGvcf
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestReblockGVCF {

    input {
      File gvcf
      File gvcf_index
      File? calling_interval_list
      File ref_dict
      File ref_fasta
      File ref_fasta_index
      Float? tree_score_cutoff
      String? annotations_to_keep_command
      String? annotations_to_remove_command
      Boolean? move_filters_to_genotypes
      String? gvcf_file_extension

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
      String cloud_provider
    }

    meta {
      allowNestedInputs: true
    }
  
    call ReblockGVCF.ReblockGVCF {
      input:
        gvcf = gvcf,
        gvcf_index = gvcf_index,
        calling_interval_list = calling_interval_list,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        tree_score_cutoff = tree_score_cutoff,
        annotations_to_keep_command = annotations_to_keep_command,
        annotations_to_remove_command = annotations_to_remove_command,
        move_filters_to_genotypes = move_filters_to_genotypes,
        gvcf_file_extension = gvcf_file_extension,
        cloud_provider = cloud_provider
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    ReblockGVCF.reblocked_gvcf_index,
                                    ReblockGVCF.reblocked_gvcf,
                                    ],
                                    
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
        call Utilities.GetValidationInputs as GetGvcf {
          input:
            input_file = ReblockGVCF.reblocked_gvcf,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetGvcfIndex {
          input:
            input_file = ReblockGVCF.reblocked_gvcf_index,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifyGvcf.VerifyGvcf as Verify {
        input:
          truth_gvcf = GetGvcf.truth_file, 
          test_gvcf = GetGvcf.results_file,
          truth_gvcf_index = GetGvcfIndex.truth_file, 
          test_gvcf_index = GetGvcfIndex.results_file,
          done = CopyToTestResults.done
      }
    }
}