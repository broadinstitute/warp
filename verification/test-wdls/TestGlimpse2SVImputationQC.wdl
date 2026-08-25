version 1.0

import "../../pipelines/wdl/glimpse/sv_imputation/input_qc/Glimpse2SVImputationQC.wdl" as Glimpse2SVImputationQC
import "../../verification/VerifyGlimpse2SVImputationQC.wdl" as VerifyGlimpse2SVImputationQC
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestGlimpse2SVImputationQC {

    input {
        File gvcf_manifest
        String output_prefix

        File preprocess_panel_bubble_split_sites_only_vcf
        File preprocess_panel_bubble_split_sites_only_vcf_idx

        Array[String] paste_regions

        Array[String] chromosomes
        File genetic_maps_tsv
        File ref_dict
        File chunked_panel_json

        File pop_glimpse2_panel_resources_json

        Float? info_filter_for_inclusion

        # for warp testing only
        String? billing_project_for_rp

        String? pipeline_header_line

        # These values will be determined and injected into the inputs by the scala test framework
        String truth_path
        String results_path
        Boolean update_truth
    }

    meta {
      allowNestedInputs: true
    }

    call Glimpse2SVImputationQC.InputQC {
      input:
        gvcf_manifest = gvcf_manifest,
        output_prefix = output_prefix,
        preprocess_panel_bubble_split_sites_only_vcf = preprocess_panel_bubble_split_sites_only_vcf,
        preprocess_panel_bubble_split_sites_only_vcf_idx = preprocess_panel_bubble_split_sites_only_vcf_idx,
        paste_regions = paste_regions,
        chromosomes = chromosomes,
        genetic_maps_tsv = genetic_maps_tsv,
        ref_dict = ref_dict,
        chunked_panel_json = chunked_panel_json,
        pop_glimpse2_panel_resources_json = pop_glimpse2_panel_resources_json,
        info_filter_for_inclusion = info_filter_for_inclusion,
        billing_project_for_rp = billing_project_for_rp,
        pipeline_header_line = pipeline_header_line
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

      call VerifyGlimpse2SVImputationQC.VerifyGlimpse2SVImputationQC as Verify {
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
