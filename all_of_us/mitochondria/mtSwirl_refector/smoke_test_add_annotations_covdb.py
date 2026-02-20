#!/usr/bin/env python3
"""Smoke test for covdb-backed mt_mean_coverage fill in add_annotations.

This test builds a tiny 2-sample MatrixTable and coverage.h5, then verifies
missing mt_mean_coverage values are filled from the coverage DB using the
same block-read strategy as finalize.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import h5py
import numpy as np
import hail as hl

from generate_mtdna_call_mt.add_annotations import (
    fill_missing_mt_mean_coverage_from_covdb,
)


def _write_covdb(path: Path) -> None:
    sample_ids = ["s1", "s2"]
    positions = np.array([1, 2, 3], dtype=np.int32)
    coverage = np.array(
        [
            [10, 20, 30],  # s1 mean=20
            [5, 5, 5],  # s2 mean=5
        ],
        dtype=np.int32,
    )

    with h5py.File(path, "w") as h5:
        str_dtype = h5py.string_dtype(encoding="utf-8")
        h5.create_dataset(
            "sample_id", data=np.array(sample_ids, dtype=object), dtype=str_dtype
        )
        h5.create_dataset("pos", data=positions)
        h5.create_dataset("coverage", data=coverage)


def _build_mt() -> hl.MatrixTable:
    mt = hl.utils.range_matrix_table(n_rows=3, n_cols=2)
    mt = mt.key_rows_by(
        locus=hl.locus("MT", mt.row_idx + 1, reference_genome="GRCh37"),
        alleles=hl.array(["A", "C"]),
    )
    mt = mt.key_cols_by(s=hl.if_else(mt.col_idx == 0, "s1", "s2"))

    mt = mt.annotate_cols(
        mt_mean_coverage=hl.if_else(mt.s == "s1", 20.0, hl.missing(hl.tfloat64))
    )
    return mt


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        covdb_path = tmp_path / "coverage.h5"
        _write_covdb(covdb_path)

        mt = _build_mt()

        mt = fill_missing_mt_mean_coverage_from_covdb(
            mt,
            coverage_h5_path=str(covdb_path),
            position_block_size=2,
        )

        results = {r.s: r.mt_mean_coverage for r in mt.cols().collect()}
        assert results["s1"] == 20.0
        assert abs(results["s2"] - 5.0) < 1e-6

        print("smoke_test_add_annotations_covdb: PASS")


if __name__ == "__main__":
    main()
