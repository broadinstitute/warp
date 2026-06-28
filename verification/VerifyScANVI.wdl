version 1.0

import "../verification/VerifyTasks.wdl" as VerifyTasks

## Verification workflow for the scANVI pipeline.
##
## scANVI trains stochastic SCVI/SCANVI models, so outputs are not bit-reproducible.
## Each output h5ad is compared to truth tolerantly (see CompareScanviH5ad): cell
## counts must match, the annotation column must be present, the predicted-label
## vocabulary must be a subset of truth's, and per-cell-type proportions must
## correlate with truth above a threshold.
##
## The ATAC-annotated output is optional: it is only produced (and only verified) in
## multiome mode. In GEX-only mode the *_atac_annotated_matrix.h5ad inputs are absent.
workflow VerifyScANVI {

    input {
        File test_scanvi_predictions_h5ad
        File truth_scanvi_predictions_h5ad

        File test_gex_annotated_h5ad
        File truth_gex_annotated_h5ad

        File? test_atac_annotated_h5ad
        File? truth_atac_annotated_h5ad

        Boolean? done
    }

    # SCANVI predictions: finalize_output renames the annotation column to 'celltype'
    call VerifyTasks.CompareScanviH5ad as CompareScanviPredictions {
        input:
            truth_h5ad = truth_scanvi_predictions_h5ad,
            test_h5ad  = test_scanvi_predictions_h5ad,
            label_key  = "celltype"
    }

    # Annotated GEX matrix: predictions stored in 'final_annotation'
    call VerifyTasks.CompareScanviH5ad as CompareGexAnnotated {
        input:
            truth_h5ad = truth_gex_annotated_h5ad,
            test_h5ad  = test_gex_annotated_h5ad,
            label_key  = "final_annotation"
    }

    # Annotated ATAC matrix: only present in multiome mode
    if (defined(test_atac_annotated_h5ad) && defined(truth_atac_annotated_h5ad)) {
        call VerifyTasks.CompareScanviH5ad as CompareAtacAnnotated {
            input:
                truth_h5ad = select_first([truth_atac_annotated_h5ad]),
                test_h5ad  = select_first([test_atac_annotated_h5ad]),
                label_key  = "final_annotation"
        }
    }

    meta {
        allowNestedInputs: true
    }
}
