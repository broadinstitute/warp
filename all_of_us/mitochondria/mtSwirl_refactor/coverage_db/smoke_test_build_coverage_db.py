#!/usr/bin/env python3
"""Smoke test for build_coverage_db.

This test is self-contained and does not require Hail.
It creates tiny synthetic coverage TSVs (two samples, five positions), runs the
builder, and asserts mean/median/over_100/over_1000 are correct.

Run manually (from repo root):
  python3 -m generate_mtdna_call_mt.coverage_db.smoke_test_build_coverage_db
"""

from __future__ import annotations

import csv
import os
import tempfile

import numpy as np

from .build_coverage_db import CoverageDBConfig, build_coverage_db, _median_int32_from_int_array


def _write_cov(path: str, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chrom", "pos", "target", "x"])
        for chrom, pos, cov in rows:
            w.writerow([chrom, pos, "NA", cov])


def _read_summary(path: str):
    out = {}
    with open(path, "r", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            pos = int(row["pos"])
            out[pos] = {
                "mean": float(row["mean"]),
                "median": int(row["median"]),
                "over_100": float(row["over_100"]),
                "over_1000": float(row["over_1000"]),
            }
    return out


def main() -> None:
    # Sanity-check median semantics against Hail docs:
    # hl.median([1, 3, 5, 6, 7, 9]) == 5 (lower-of-two-middles for even N)
    assert int(_median_int32_from_int_array(np.array([1, 3, 5, 6, 7, 9], dtype=np.int64))) == 5

    with tempfile.TemporaryDirectory() as d:
        cov1 = os.path.join(d, "s1.tsv")
        cov2 = os.path.join(d, "s2.tsv")

        # positions 1..5
        # s1: [0, 50, 101, 1001, 7]
        # s2: [0, 150, 99,  999,  7]
        _write_cov(cov1, [("chrM", 1, 0), ("chrM", 2, 50), ("chrM", 3, 101), ("chrM", 4, 1001), ("chrM", 5, 7)])
        _write_cov(cov2, [("chrM", 1, 0), ("chrM", 2, 150), ("chrM", 3, 99), ("chrM", 4, 999), ("chrM", 5, 7)])

        out_h5 = os.path.join(d, "coverage.h5")
        out_sum = os.path.join(d, "coverage_summary.tsv")

        cfg = CoverageDBConfig(position_start=1, position_end=5, chrom="chrM")

        build_coverage_db(
            samples=[("s1", cov1), ("s2", cov2)],
            out_h5=out_h5,
            out_summary_tsv=out_sum,
            config=cfg,
            batch_size=2,
            position_block_size=3,
            compression=None,
            compression_level=None,
            dtype="uint16",
            expected_chrom="chrM",
            skip_summary=False,
        )

        got = _read_summary(out_sum)

        # Expected:
        # pos1: mean=0 median=0 over_100=0 over_1000=0
    # pos2: mean=100 median=50 over_100=0.5 (150>100) over_1000=0
    # pos3: mean=100 median=99 over_100=0.5 (101>100) over_1000=0
    # pos4: mean=1000 median=999 over_100=1 over_1000=0.5 (1001>1000)
        # pos5: mean=7 median=7 over_100=0 over_1000=0
        exp = {
            1: {"mean": 0.0, "median": 0, "over_100": 0.0, "over_1000": 0.0},
            2: {"mean": 100.0, "median": 50, "over_100": 0.5, "over_1000": 0.0},
            3: {"mean": 100.0, "median": 99, "over_100": 0.5, "over_1000": 0.0},
            4: {"mean": 1000.0, "median": 999, "over_100": 1.0, "over_1000": 0.5},
            5: {"mean": 7.0, "median": 7, "over_100": 0.0, "over_1000": 0.0},
        }

        for pos, e in exp.items():
            g = got[pos]
            assert np.isclose(g["mean"], e["mean"])
            assert g["median"] == e["median"]
            assert np.isclose(g["over_100"], e["over_100"])
            assert np.isclose(g["over_1000"], e["over_1000"])

        print("OK: smoke test passed")


if __name__ == "__main__":
    main()
