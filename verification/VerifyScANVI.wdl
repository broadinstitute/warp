version 1.0

import "../tasks/wdl/Utilities.wdl" as Utilities

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
    call CompareScanviH5ad as CompareScanviPredictions {
        input:
            truth_h5ad = truth_scanvi_predictions_h5ad,
            test_h5ad  = test_scanvi_predictions_h5ad,
            label_key  = "celltype"
    }

    # Annotated GEX matrix: predictions stored in 'final_annotation'
    call CompareScanviH5ad as CompareGexAnnotated {
        input:
            truth_h5ad = truth_gex_annotated_h5ad,
            test_h5ad  = test_gex_annotated_h5ad,
            label_key  = "final_annotation"
    }

    # ATAC-annotated output is optional (produced only in multiome mode). Fail loudly if
    # exactly one side defines it — that signals an unexpected presence/absence regression
    # (e.g. truth has ATAC but a code change made the run GEX-only) that would otherwise be
    # silently skipped.
    if (defined(test_atac_annotated_h5ad) != defined(truth_atac_annotated_h5ad)) {
        call Utilities.ErrorWithMessage as AtacPresenceMismatch {
            input:
                message = "scANVI ATAC-annotated output is present on only one of {test, truth}; expected both (multiome) or neither (GEX-only)."
        }
    }

    # Compare ATAC when both sides have it (multiome).
    if (defined(test_atac_annotated_h5ad) && defined(truth_atac_annotated_h5ad)) {
        call CompareScanviH5ad as CompareAtacAnnotated {
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

# Tolerant comparison for the scANVI pipeline outputs. SCVI/SCANVI training is
# stochastic, so predicted labels/embeddings/counts are NOT bit-reproducible
# across runs. Instead of exact equality, this verifies that test and truth are
# structurally and distributionally consistent:
#   - same number of cells (n_obs)
#   - the annotation column (label_key) is present in both
#   - the test predicted-label vocabulary is a subset of truth's (no novel labels)
#   - per-cell-type proportions correlate with truth at or above min_proportion_corr
#
# Kept in this scANVI-specific file (not the shared verification/VerifyTasks.wdl) so
# that scANVI changes do not trigger every pipeline's path-filtered CI.
task CompareScanviH5ad {
  input {
    File truth_h5ad
    File test_h5ad
    String label_key
    Float min_proportion_corr = 0.95
    String docker = "python:3.10.0-buster"
    Int disk_size_gb = ceil(size(truth_h5ad, "GiB") + size(test_h5ad, "GiB")) + 50
    Int memory_gb = 16
  }

  command <<<
    set -eo pipefail

    pip3 install --quiet anndata

    python3 <<CODE
    import sys
    import anndata as ad
    import numpy as np
    import pandas as pd

    label_key = "~{label_key}"
    min_corr = float("~{min_proportion_corr}")

    truth = ad.read_h5ad("~{truth_h5ad}")
    test = ad.read_h5ad("~{test_h5ad}")
    print(f"truth: {truth.shape} | test: {test.shape} | label_key={label_key}")

    # 1. Cell counts must match truth
    if truth.n_obs != test.n_obs:
        sys.exit(f"FAIL: n_obs differs (truth={truth.n_obs}, test={test.n_obs})")

    # 2. Annotation column must be present in both
    for name, a in [("truth", truth), ("test", test)]:
        if label_key not in a.obs.columns:
            sys.exit(f"FAIL: label_key '{label_key}' missing from {name}.obs")

    truth_labels = truth.obs[label_key].astype(str)
    test_labels = test.obs[label_key].astype(str)

    # 3. No novel labels in test that the reference (truth) never produced
    novel = set(test_labels.unique()) - set(truth_labels.unique())
    if novel:
        sys.exit(f"FAIL: test produced labels not present in truth: {sorted(novel)}")
    else:
        print("PASS: test label vocabulary is a subset of truth's")

    # 4. Per-cell-type proportions must correlate with truth
    all_labels = sorted(set(truth_labels.unique()) | set(test_labels.unique()))
    tp = truth_labels.value_counts(normalize=True).reindex(all_labels).fillna(0.0).to_numpy()
    sp = test_labels.value_counts(normalize=True).reindex(all_labels).fillna(0.0).to_numpy()
    if len(all_labels) < 2:
        corr = 1.0  # degenerate single-label case
    else:
        with np.errstate(invalid="ignore"):
            corr = float(np.corrcoef(tp, sp)[0, 1])
        # corrcoef is NaN when a proportion vector has zero variance (e.g. uniform
        # distribution): treat identical distributions as a perfect match, otherwise
        # as no correlation so the threshold check can fail rather than pass on NaN.
        if np.isnan(corr):
            corr = 1.0 if np.allclose(tp, sp) else 0.0
    print("Per-cell-type proportions (truth vs test):")
    for lab, t, s in zip(all_labels, tp, sp):
        print(f"  {lab:40s} truth={t:.4f}  test={s:.4f}")
    print(f"proportion correlation = {corr:.4f} (threshold {min_corr})")
    if corr < min_corr:
        sys.exit(f"FAIL: cell-type proportion correlation {corr:.4f} < {min_corr}")

    print("PASS: scANVI outputs are structurally and distributionally consistent")
    CODE
  >>>

  runtime {
    docker: docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_gb} GiB"
    preemptible: 3
  }
}
