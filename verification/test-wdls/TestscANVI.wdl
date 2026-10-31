version 1.0


import "../../pipelines/wdl/scanvi/scANVI.wdl" as scANVI
import "../../verification/VerifyScANVI.wdl" as VerifyScANVI
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestscANVI {

    input {
      # scANVI pipeline inputs: query arrays paired by position (a single query is a one-element
      # array). ATAC is optional; leave atac_h5ads empty for GEX-only mode.
      Array[String] input_ids
      Array[File] gex_h5ads
      File? ref_h5ad
      Array[File] atac_h5ads = []

      # Optional reference-column / genome overrides, forwarded to scANVI. Unset => the
      # pipeline's own defaults (AIT subclass/donor_id; genome hg38).
      Array[String] ref_label_columns = []
      String? ref_batch_column
      String genome = "hg38"
      Boolean output_max_probability = true

      # Optional cap on SCVI/SCANVI training epochs, forwarded to scANVI (the training-Plumbing
      # test sets this to 2 for a fast run; unset => the pipeline's container default of 500).
      Int? max_epochs
      # Epochs when mapping a sample onto a reference-only model, forwarded to scANVI.
      Int? query_max_epochs

      # SCVI/SCANVI minibatch size, forwarded to scANVI. Default 128. Exposed for completeness /
      # large-VRAM cloud runs; our tests leave it at the default.
      Int batch_size = 128
      # Fixed seed so SCVI/SCANVI training is reproducible against truth.
      Int seed = 42

      # Optional pre-trained SCANVI model (skip training) + compute overrides, forwarded to scANVI.
      # The pretrained Plumbing test sets scanvi_model + gpu_count=0 + small mem/cpu/disk.
      File? scanvi_model
      Int gpu_count = 2
      Int mem_size = 104
      Int nthreads = 16
      Int? disk_size

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
      # Skip the label-distribution checks only for Plumbing runs that train (a 2-epoch model is
      # meant to vary); Scientific runs and any run that loads a fixed scanvi_model keep them.
      # Plumbing vs Scientific comes from the framework's truth_path (.../truth/<plumbing|scientific>/...).
      Boolean check_label_distribution = defined(scanvi_model) || sub(truth_path, "/truth/plumbing/", "") == truth_path
    }

    meta {
      allowNestedInputs: true
    }

    call scANVI.scANVI {
      input:
        input_ids        = input_ids,
        gex_h5ads        = gex_h5ads,
        ref_h5ad         = ref_h5ad,
        atac_h5ads       = atac_h5ads,
        ref_label_columns = ref_label_columns,
        ref_batch_column = ref_batch_column,
        genome           = genome,
        output_max_probability = output_max_probability,
        max_epochs       = max_epochs,
        query_max_epochs = query_max_epochs,
        batch_size       = batch_size,
        seed             = seed,
        scanvi_model     = scanvi_model,
        gpu_count        = gpu_count,
        mem_size         = mem_size,
        nthreads         = nthreads,
        disk_size        = disk_size
    }

    # Collect all of the pipeline outputs (every query's) into a single Array[String].
    # atac_annotated_h5ads is empty in GEX-only mode, scanvi_models_out when a model is supplied.
    # scanvi_models_out is archived to results/truth (not verified — model weights are
    # non-deterministic across training runs).
    Array[String] pipeline_outputs = flatten([
                                    scANVI.gex_annotated_h5ads,
                                    scANVI.scanvi_models_out,
                                    scANVI.atac_annotated_h5ads,
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
      # Verify each query's outputs against its own truth files (matched by file name; outputs
      # are named after the input files, which must have distinct names).
      scatter (i in range(length(input_ids))) {
        call Utilities.GetValidationInputs as GetGexAnnotated {
          input:
            input_file = scANVI.gex_annotated_h5ads[i],
            results_path = results_path,
            truth_path = truth_path
        }
        # ATAC-annotated output only exists in multiome mode (keyed off the test's own input, so a
        # multiome test that loses its ATAC outputs fails here instead of skipping the comparison)
        if (length(atac_h5ads) > 0) {
          call Utilities.GetValidationInputs as GetAtacAnnotated {
            input:
              input_file = scANVI.atac_annotated_h5ads[i],
              results_path = results_path,
              truth_path = truth_path
          }
        }

        call VerifyScANVI.VerifyScANVI as Verify {
          input:
            truth_gex_annotated_h5ad      = GetGexAnnotated.truth_file,
            test_gex_annotated_h5ad       = GetGexAnnotated.results_file,
            truth_atac_annotated_h5ad     = GetAtacAnnotated.truth_file,
            test_atac_annotated_h5ad      = GetAtacAnnotated.results_file,
            check_label_distribution      = check_label_distribution,
            expected_input_id             = input_ids[i],
            done = CopyToTestResults.done
        }
      }
    }
}
