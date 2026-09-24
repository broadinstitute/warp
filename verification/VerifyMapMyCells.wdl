version 1.0

## Verification workflow for the MapMyCells pipeline.
##
## cell_type_mapper's bootstrapping step resamples marker genes, so per-cell
## probability columns are not guaranteed bit-reproducible across runs even with
## identical inputs. CompareMapMyCellsOutputs therefore checks structural /
## distributional agreement rather than exact equality: same cell set, matching
## hierarchical label assignments above an agreement threshold, and bootstrap
## probabilities within tolerance. Column names are not assumed beyond a leading
## cell-id column and `*_label` / `*_bootstrapping_probability` suffixes (the
## exact schema hasn't been pinned against a real run yet -- see MapMyCells
## pipeline README/changelog for the outstanding reference-data blocker).
workflow VerifyMapMyCells {

    input {
        File test_output_csv
        File truth_output_csv

        File test_output_json
        File truth_output_json

        Boolean? done
    }

    call CompareMapMyCellsOutputs {
        input:
            truth_csv  = truth_output_csv,
            test_csv   = test_output_csv,
            truth_json = truth_output_json,
            test_json  = test_output_json
    }

    meta {
        allowNestedInputs: true
    }
}

task CompareMapMyCellsOutputs {
  input {
    File truth_csv
    File test_csv
    File truth_json
    File test_json
    Float min_label_agreement = 0.95
    Float max_mean_probability_diff = 0.1
    String docker = "python:3.10.0-buster"
    Int disk_size_gb = ceil(size(truth_csv, "GiB") + size(test_csv, "GiB") + size(truth_json, "GiB") + size(test_json, "GiB")) + 10
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    pip3 install --quiet pandas

    python3 <<CODE
    import json
    import sys
    import pandas as pd

    min_agreement = float("~{min_label_agreement}")
    max_prob_diff = float("~{max_mean_probability_diff}")

    truth = pd.read_csv("~{truth_csv}", comment="#")
    test = pd.read_csv("~{test_csv}", comment="#")
    print(f"truth csv: {truth.shape} | test csv: {test.shape}")

    if truth.shape[0] != test.shape[0]:
        sys.exit(f"FAIL: row count differs (truth={truth.shape[0]}, test={test.shape[0]})")

    id_col = truth.columns[0]
    if id_col != test.columns[0]:
        sys.exit(f"FAIL: leading (cell id) column differs (truth={id_col}, test={test.columns[0]})")

    truth = truth.set_index(id_col).sort_index()
    test = test.set_index(id_col).sort_index()
    if not truth.index.equals(test.index):
        sys.exit("FAIL: cell id sets differ between truth and test")

    label_cols = [c for c in truth.columns if c.endswith("_label")]
    if not label_cols:
        sys.exit("FAIL: no '*_label' columns found in truth csv")

    for col in label_cols:
        if col not in test.columns:
            sys.exit(f"FAIL: label column '{col}' missing from test csv")
        agreement = (truth[col].astype(str) == test[col].astype(str)).mean()
        print(f"{col}: agreement={agreement:.4f} (threshold {min_agreement})")
        if agreement < min_agreement:
            sys.exit(f"FAIL: '{col}' agreement {agreement:.4f} < {min_agreement}")

    prob_cols = [c for c in truth.columns if c.endswith("_bootstrapping_probability")]
    for col in prob_cols:
        if col not in test.columns:
            sys.exit(f"FAIL: probability column '{col}' missing from test csv")
        mean_diff = (truth[col] - test[col]).abs().mean()
        print(f"{col}: mean abs diff={mean_diff:.4f} (threshold {max_prob_diff})")
        if mean_diff > max_prob_diff:
            sys.exit(f"FAIL: '{col}' mean abs diff {mean_diff:.4f} > {max_prob_diff}")

    with open("~{truth_json}") as f:
        truth_json = json.load(f)
    with open("~{test_json}") as f:
        test_json = json.load(f)

    def key_shape(obj):
        if isinstance(obj, dict):
            return {k: key_shape(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [key_shape(obj[0])] if obj else []
        return None

    if key_shape(truth_json) != key_shape(test_json):
        sys.exit("FAIL: extended-result JSON key structure differs between truth and test")

    print("PASS: MapMyCells outputs are structurally and distributionally consistent")
    CODE
  >>>

  runtime {
    docker: docker
    disks: "local-disk ~{disk_size_gb} HDD"
    memory: "~{memory_gb} GiB"
    preemptible: 3
  }
}
