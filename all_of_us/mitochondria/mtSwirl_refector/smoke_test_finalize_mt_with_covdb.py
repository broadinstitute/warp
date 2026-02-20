#!/usr/bin/env python3
"""Local smoke test for covdb-based hom-ref imputation.

This is a tiny, deterministic test that creates:
- a 2-sample, 3-position Hail MatrixTable
- a matching `coverage.h5`

It then runs `determine_hom_refs_from_covdb` and asserts that DP/HL/FT
match the expected v1 semantics.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import h5py
import numpy as np
import hail as hl

from generate_mtdna_call_mt.merging_utils import (
    apply_mito_artifact_filter,
    determine_hom_refs_from_covdb,
)


def _write_covdb(path: Path) -> None:
    sample_ids = ["s1", "s2"]
    positions = np.array([1, 2, 3], dtype=np.int32)
    coverage = np.array(
        [
            [200, 50, 150],  # s1
            [10, 200, 0],  # s2
        ],
        dtype=np.int32,
    )

    with h5py.File(path, "w") as h5:
        str_dtype = h5py.string_dtype(encoding="utf-8")
        h5.create_dataset("sample_id", data=np.array(sample_ids, dtype=object), dtype=str_dtype)
        h5.create_dataset("pos", data=positions)
        h5.create_dataset("coverage", data=coverage)


def _build_mt() -> hl.MatrixTable:
    mt = hl.utils.range_matrix_table(n_rows=3, n_cols=2)
    mt = mt.key_rows_by(
        locus=hl.locus("MT", mt.row_idx + 1, reference_genome="GRCh37"),
        alleles=hl.array(["A", "C"]),
    )
    mt = mt.key_cols_by(s=hl.if_else(mt.col_idx == 0, "s1", "s2"))

    keep_entry = (mt.row_idx == 1) & (mt.col_idx == 0)
    mt = mt.annotate_entries(
        HL=hl.or_missing(keep_entry, 0.5),
        DP=hl.or_missing(keep_entry, 20),
        FT=hl.set(["PASS"]),
    )
    return mt


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        covdb_path = tmp_path / "coverage.h5"
        _write_covdb(covdb_path)

        bed_path = tmp_path / "artifact.bed"
        bed_path.write_text("chrM\t1\t2\n")

        hl.init(tmp_dir=str(tmp_path / "hail_tmp"), quiet=True)
        mt = _build_mt()

        mt = determine_hom_refs_from_covdb(
            mt,
            coverage_h5_path=str(covdb_path),
            minimum_homref_coverage=100,
            position_block_size=2,
        )

        mt = apply_mito_artifact_filter(
            mt,
            artifact_prone_sites_path=str(bed_path),
            artifact_prone_sites_reference="default",
        )

        mt = mt.checkpoint(str(tmp_path / "final.mt"), overwrite=True, stage_locally=True)

        ht_entries = mt.entries()
        entries = (
            ht_entries.select(
                pos=ht_entries.locus.position,
                HL=ht_entries.HL,
                DP=ht_entries.DP,
                FT=ht_entries.FT,
            )
            .collect()
        )

        result = {(e.s, int(e.pos)): (e.HL, e.DP) for e in entries}

        # s1
        assert result[("s1", 1)] == (0.0, 200)
        assert result[("s1", 2)] == (0.5, 20)
        assert result[("s1", 3)] == (0.0, 150)
        # s2
        assert result[("s2", 1)] == (None, None)
        assert result[("s2", 2)] == (0.0, 200)
        assert result[("s2", 3)] == (None, None)

        print("smoke_test_finalize_mt_with_covdb: PASS")


if __name__ == "__main__":
    main()
