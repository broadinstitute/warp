"""Smoke test for covdb-based hom-ref / DP imputation.

This test focuses on the core contract replicated from v1 `determine_hom_refs`.
The implementation under test is v2 and uses a block-broadcast strategy internally.
- For missing HL, pull DP from coverage.
- If DP > minimum_homref_coverage => HL=0.0 and FT={"PASS"}.
- If DP <= minimum_homref_coverage => DP becomes missing, HL remains missing.

We avoid full multi-way MT merging and heavy dependencies; instead we create a tiny
MatrixTable in Hail directly and apply `determine_hom_refs_from_covdb`.

Run (in an environment with hail + h5py):
    python -m generate_mtdna_call_mt.smoke_test_combine_vcfs_and_homref_from_covdb
"""

from __future__ import annotations

import os
import tempfile

import hail as hl
import numpy as np

try:
    import h5py  # type: ignore
except Exception as e:
    raise RuntimeError("h5py is required for this smoke test") from e

from generate_mtdna_call_mt.merging_utils import determine_hom_refs_from_covdb


def _write_tiny_covdb(h5_path: str) -> None:
    # Two samples, three positions
    sample_id = np.array(["S1", "S2"], dtype=object)
    pos = np.array([1, 2, 3], dtype=np.int32)
    # coverage shape: (n_samples, n_positions)
    # S1: high coverage at pos1, low at pos2, high at pos3
    # S2: low at pos1, high at pos2, low at pos3
    coverage = np.array(
        [
            [200, 50, 150],
            [10, 300, 20],
        ],
        dtype=np.uint16,
    )

    with h5py.File(h5_path, "w") as h5:
        h5.create_dataset("coverage", data=coverage, compression=None)
        # vlen string dataset
        dt = h5py.string_dtype(encoding="utf-8")
        h5.create_dataset("sample_id", data=sample_id.astype(dt), dtype=dt)
        h5.create_dataset("pos", data=pos)


def main() -> None:
    hl.init(quiet=True)

    with tempfile.TemporaryDirectory() as tmp:
        h5_path = os.path.join(tmp, "coverage.h5")
        _write_tiny_covdb(h5_path)

        # Build a tiny MT with 2 samples x 3 rows and HL missing everywhere.
        mt = hl.utils.range_matrix_table(n_rows=3, n_cols=2)
        mt = mt.annotate_rows(locus=hl.locus("MT", mt.row_idx + 1, reference_genome="GRCh37"), alleles=["A", "C"]).key_rows_by("locus", "alleles")
        mt = mt.annotate_cols(s=hl.if_else(mt.col_idx == 0, "S1", "S2")).key_cols_by("s")

        # Entries: HL missing, DP present but should be overwritten only when HL missing (always true).
        mt = mt.annotate_entries(HL=hl.missing(hl.tfloat64), DP=hl.missing(hl.tint32), FT={"FAIL"})

        mt2 = determine_hom_refs_from_covdb(
            mt,
            coverage_h5_path=h5_path,
            minimum_homref_coverage=100,
            position_block_size=2,
        )

        # Collect results into python for assertions.
        ht = mt2.entries().key_by("locus", "s")
        res = {(r.locus.position, r.s): (r.DP, r.HL, set(r.FT)) for r in ht.collect()}

        # Position 1: S1 cov=200 => HL=0, DP=200, FT=PASS; S2 cov=10 => DP missing, HL missing
        assert res[(1, "S1")][0] == 200
        assert res[(1, "S1")][1] == 0.0
        assert res[(1, "S1")][2] == {"PASS"}
        assert res[(1, "S2")][0] is None
        assert res[(1, "S2")][1] is None

        # Position 2: S1 cov=50 => DP missing, HL missing; S2 cov=300 => HL=0, DP=300
        assert res[(2, "S1")][0] is None
        assert res[(2, "S1")][1] is None
        assert res[(2, "S2")][0] == 300
        assert res[(2, "S2")][1] == 0.0
        assert res[(2, "S2")][2] == {"PASS"}

        # Position 3: S1 cov=150 => HL=0; S2 cov=20 => missing
        assert res[(3, "S1")][0] == 150
        assert res[(3, "S1")][1] == 0.0
        assert res[(3, "S2")][0] is None
        assert res[(3, "S2")][1] is None

        print("OK: covdb hom-ref smoke test passed")


if __name__ == "__main__":
    main()
