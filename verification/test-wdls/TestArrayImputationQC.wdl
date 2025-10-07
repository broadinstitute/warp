version 1.0


import "../../pipelines/wdl/arrays/imputation_beagle/ArrayImputationQC.wdl" as ArrayImputationQC
import "../../verification/VerifyArrayImputationQC.wdl" as VerifyArrayImputationQC
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestArrayImputationQC {

    input {
      Int chunkLength = 25000000
      Int chunkOverlaps = 5000000 # this is the padding that will be added to the beginning and end of each chunk to reduce edge effects
      Int sample_chunk_size = 1000 # this is the number of samples that will be processed in parallel in each chunked scatter
      Float min_dr2_for_inclusion = 0.0 # minimum imputation quality (DR2) for a variant to be included in the output VCF
      
      File multi_sample_vcf
      
      File ref_dict # for reheadering / adding contig lengths in the header of the ouptut VCF, and calculating contig lengths
      Array[String] contigs # list of possible contigs that will be processed. note the workflow will not error out if any of these contigs are missing
      String reference_panel_path_prefix # path + file prefix to the bucket where the reference panel files are stored for all contigs
      String genetic_maps_path # path to the bucket where genetic maps are stored for all contigs
      String output_basename # the basename for intermediate and output files

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
    }

    meta {
      allowNestedInputs: true
    }
  
    call ArrayImputationQC.InputQC {
      input:
        chunkLength = chunkLength,
        chunkOverlaps = chunkOverlaps,
        multi_sample_vcf = multi_sample_vcf,
        ref_dict = ref_dict,
        contigs = contigs,
        reference_panel_path_prefix = reference_panel_path_prefix,
        genetic_maps_path = genetic_maps_path,
        output_basename = output_basename,
        min_dr2_for_inclusion = min_dr2_for_inclusion,
    }

    # Write pipeline outputs into json file so we can compare to truth
    call WritePipelineOutputs {
      input:
        input_map = {
          "passes_qc": InputQC.passes_qc,
          "qc_messages": InputQC.qc_messages
        }
    }
    
    # Copy results of pipeline to test results bucket
    call Copy.TerraCopyFilesFromCloudToCloud as CopyToTestResults {
      input:
        files_to_copy             = [WritePipelineOutputs.json_file],
        destination_cloud_path    = results_path
    }
  
    # If updating truth then copy output to truth bucket
    if (update_truth){
      call Copy.TerraCopyFilesFromCloudToCloud as CopyToTruth {
        input: 
          files_to_copy             = [WritePipelineOutputs.json_file],
          destination_cloud_path    = truth_path
      }
    }

    # This is achieved by passing each desired file/array[files] to GetValidationInputs
    if (!update_truth){
        call Utilities.GetValidationInputs as GetOutputs {
          input:
            input_file = WritePipelineOutputs.json_file,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifyArrayImputationQC.VerifyArrayImputationQC as Verify {
        input:
          truth_outputs = GetOutputs.truth_file, 
          test_outputs = GetOutputs.results_file,
          done = CopyToTestResults.done
      }
    }
}

# Write a json file from a map of strings
task WritePipelineOutputs {
  input {
    Map[String, String] input_map
  }

  command <<<
    python3 -c 'import json; print(json.dumps(~{input_map}, indent=4))' > output.json
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
    memory: "2 GiB"
  }

  output {
    File json_file = "output.json"
  }
}
