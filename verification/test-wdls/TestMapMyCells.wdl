version 1.0

import "../../pipelines/wdl/mapmycells/MapMyCells.wdl" as MapMyCells
import "../../verification/VerifyMapMyCells.wdl" as VerifyMapMyCells
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestMapMyCells {

    input {
      # MapMyCells pipeline inputs
      File query_h5ad
      String input_id
      String reference_atlas = "Human_MTG"
      File? custom_precomputed_stats
      File? custom_gene_mapping_db
      File? custom_query_markers
      String algorithm = "hierarchical"
      String normalization = "raw"
      Int cpu = 8
      Int memory_gb = 64
      Int disk_size = 100
      Int preemptible = 3

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
    }

    meta {
      allowNestedInputs: true
    }

    call MapMyCells.MapMyCells {
      input:
        query_h5ad               = query_h5ad,
        input_id                 = input_id,
        reference_atlas          = reference_atlas,
        custom_precomputed_stats = custom_precomputed_stats,
        custom_gene_mapping_db   = custom_gene_mapping_db,
        custom_query_markers     = custom_query_markers,
        algorithm                = algorithm,
        normalization            = normalization,
        cpu                      = cpu,
        memory_gb                = memory_gb,
        disk_size                = disk_size,
        preemptible               = preemptible
    }

    Array[String] pipeline_outputs = [
        MapMyCells.output_csv,
        MapMyCells.output_json,
    ]

    # Copy results of pipeline to test results bucket
    call Copy.TerraCopyFilesFromCloudToCloud as CopyToTestResults {
      input:
        files_to_copy             = pipeline_outputs,
        destination_cloud_path    = results_path
    }

    # If updating truth then copy output to truth bucket
    if (update_truth){
      call Copy.TerraCopyFilesFromCloudToCloud as CopyToTruth {
        input:
          files_to_copy             = pipeline_outputs,
          destination_cloud_path    = truth_path
      }
    }

    if (!update_truth){
        call Utilities.GetValidationInputs as GetOutputCsv {
          input:
            input_file = MapMyCells.output_csv,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetOutputJson {
          input:
            input_file = MapMyCells.output_json,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifyMapMyCells.VerifyMapMyCells as Verify {
        input:
          truth_output_csv  = GetOutputCsv.truth_file,
          test_output_csv   = GetOutputCsv.results_file,
          truth_output_json = GetOutputJson.truth_file,
          test_output_json  = GetOutputJson.results_file,
          done = CopyToTestResults.done
      }
    }
}
