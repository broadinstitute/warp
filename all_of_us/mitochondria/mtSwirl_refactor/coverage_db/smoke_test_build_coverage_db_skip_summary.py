#!/usr/bin/env python3
"""Smoke test for build_coverage_db with --skip-summary.

This verifies that when summary is skipped we still produce a valid coverage.h5
and we do NOT create coverage_summary.tsv.

Run manually (from repo root):
  python3 -m generate_mtdna_call_mt.coverage_db.smoke_test_build_coverage_db_skip_summary
"""

from __future__ import annotations

import csv
import os
import tempfile

from .build_coverage_db import CoverageDBConfig, build_coverage_db


def _write_cov(path: str, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chrom", "pos", "target", "x"])
        for chrom, pos, cov in rows:
            w.writerow([chrom, pos, "NA", cov])


def main() -> None:
    with tempfile.TemporaryDirectory() as d:
        cov1 = os.path.join(d, "s1.tsv")
        cov2 = os.path.join(d, "s2.tsv")

        _write_cov(cov1, [("chrM", 1, 0), ("chrM", 2, 50), ("chrM", 3, 101)])
        _write_cov(cov2, [("chrM", 1, 0), ("chrM", 2, 150), ("chrM", 3, 99)])

        out_h5 = os.path.join(d, "coverage.h5")
        out_sum = os.path.join(d, "coverage_summary.tsv")

        cfg = CoverageDBConfig(position_start=1, position_end=3, chrom="chrM")

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
            skip_summary=True,
        )

        assert os.path.exists(out_h5), "coverage.h5 was not created"
        assert not os.path.exists(out_sum), "coverage_summary.tsv should not be created when skip_summary=True"

        print("OK: skip_summary smoke test passed")


if __name__ == "__main__":
    main()
