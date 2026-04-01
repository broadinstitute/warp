version 1.0

import "../../pipelines/wdl/glimpse/low_pass_imputation/input_qc/LowPassImputationQC.wdl" as LowPassImputationQC
import "../../verification/VerifyLowPassImputationQC.wdl" as VerifyLowPassImputationQC
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestLowPassImputationQC {

    input {
        String reference_panel_prefix

        Array[String] contigs

        Array[File]? crams
        Array[File]? cram_indices
        Array[String]? sample_ids
        File? cram_manifest

        File fasta
        File fasta_index
        String output_basename

        File ref_dict

        # These values will be determined and injected into the inputs by the scala test framework
        String truth_path
        String results_path
        Boolean update_truth
    }

    meta {
      allowNestedInputs: true
    }
  
    call LowPassImputationQC.InputQC {
      input:
        contigs = contigs,
        reference_panel_prefix = reference_panel_prefix,
        crams = crams,
        cram_indices = cram_indices,
        sample_ids = sample_ids,
        cram_manifest = cram_manifest,
        fasta = fasta,
        fasta_index = fasta_index,
        output_basename = output_basename,
        ref_dict = ref_dict
    }

    # Write pipeline outputs into json file so we can compare to truth
    call WriteMapToTsv {
      input:
        input_map = {
          "passes_qc": InputQC.passes_qc,
          "qc_messages": InputQC.qc_messages
        }
    }
    
    # Copy results of pipeline to test results bucket
    call Copy.TerraCopyFilesFromCloudToCloud as CopyToTestResults {
      input:
        files_to_copy             = [WriteMapToTsv.tsv_file],
        destination_cloud_path    = results_path
    }
  
    # If updating truth then copy output to truth bucket
    if (update_truth){
      call Copy.TerraCopyFilesFromCloudToCloud as CopyToTruth {
        input: 
          files_to_copy             = [WriteMapToTsv.tsv_file],
          destination_cloud_path    = truth_path
      }
    }

    # This is achieved by passing each desired file/array[files] to GetValidationInputs
    if (!update_truth){
        call Utilities.GetValidationInputs as GetOutputs {
          input:
            input_file = WriteMapToTsv.tsv_file,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifyLowPassImputationQC.VerifyLowPassImputationQC as Verify {
        input:
          truth_outputs = GetOutputs.truth_file, 
          test_outputs = GetOutputs.results_file,
          done = CopyToTestResults.done
      }
    }
}

# Write a tsv file from a map of strings
task WriteMapToTsv {
  input {
    Map[String, String] input_map
  }

  command <<<
    cp ~{write_map(input_map)} output.tsv
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
    memory: "2 GiB"
  }

  output {
    File tsv_file = "output.tsv"
  }
}
