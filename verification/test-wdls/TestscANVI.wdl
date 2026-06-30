version 1.0


import "../../pipelines/wdl/scanvi/scANVI.wdl" as scANVI
import "../../verification/VerifyScANVI.wdl" as VerifyScANVI
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestscANVI {

    input {
      # scANVI pipeline inputs (ATAC is optional; omit it for GEX-only mode)
      String? input_bucket
      File? gex_h5ad
      File? atac_h5ad
      File? ref_h5ad
      String gex_filename = "gex.h5ad"
      String atac_filename = "atac.h5ad"
      String ref_filename = "ref.h5ad"
      String input_id

      # Optional reference-column / genome overrides, forwarded to scANVI. Unset => the
      # pipeline's own defaults (AIT subclass/donor_id; genome hg38), so existing tests are unchanged.
      String? ref_label_column
      String? ref_batch_column
      String genome = "hg38"
      Boolean output_max_probability = false

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
    }

    meta {
      allowNestedInputs: true
    }

    call scANVI.scANVI {
      input:
        input_bucket  = input_bucket,
        gex_h5ad      = gex_h5ad,
        atac_h5ad     = atac_h5ad,
        ref_h5ad      = ref_h5ad,
        gex_filename  = gex_filename,
        atac_filename = atac_filename,
        ref_filename     = ref_filename,
        input_id         = input_id,
        ref_label_column = ref_label_column,
        ref_batch_column = ref_batch_column,
        genome           = genome,
        output_max_probability = output_max_probability
    }

    # Collect all of the pipeline outputs into a single Array[String].
    # atac_annotated_h5ad is optional (only produced in multiome mode), so it is
    # included via select_all and contributes nothing in GEX-only mode.
    Array[String] pipeline_outputs = flatten([
                                    [ # always produced
                                    scANVI.scanvi_predictions_h5ad,
                                    scANVI.gex_annotated_h5ad,
                                    ],
                                    select_all([scANVI.atac_annotated_h5ad]),
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

    if (!update_truth){
        call Utilities.GetValidationInputs as GetScanviPredictions {
          input:
            input_file = scANVI.scanvi_predictions_h5ad,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetGexAnnotated {
          input:
            input_file = scANVI.gex_annotated_h5ad,
            results_path = results_path,
            truth_path = truth_path
        }
        # ATAC-annotated output only exists in multiome mode
        if (defined(scANVI.atac_annotated_h5ad)) {
          call Utilities.GetValidationInputs as GetAtacAnnotated {
            input:
              input_file = select_first([scANVI.atac_annotated_h5ad]),
              results_path = results_path,
              truth_path = truth_path
          }
        }

      call VerifyScANVI.VerifyScANVI as Verify {
        input:
          truth_scanvi_predictions_h5ad = GetScanviPredictions.truth_file,
          test_scanvi_predictions_h5ad  = GetScanviPredictions.results_file,
          truth_gex_annotated_h5ad      = GetGexAnnotated.truth_file,
          test_gex_annotated_h5ad       = GetGexAnnotated.results_file,
          truth_atac_annotated_h5ad     = GetAtacAnnotated.truth_file,
          test_atac_annotated_h5ad      = GetAtacAnnotated.results_file,
          done = CopyToTestResults.done
      }
    }
}
