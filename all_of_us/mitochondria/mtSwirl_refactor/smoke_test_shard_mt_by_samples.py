#!/usr/bin/env python3
"""Smoke test for shard_mt_by_samples.

Creates a small MT, shards it by samples, and validates shard sizes.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import hail as hl

import argparse

from generate_mtdna_call_mt.Terra.shard_mt_by_samples import main as shard_main


def _build_mt() -> hl.MatrixTable:
    mt = hl.utils.range_matrix_table(n_rows=2, n_cols=5)
    mt = mt.key_rows_by(
        locus=hl.locus("MT", mt.row_idx + 1, reference_genome="GRCh37"),
        alleles=hl.array(["A", "C"]),
    )
    mt = mt.key_cols_by(s=hl.str(mt.col_idx))
    mt = mt.annotate_entries(HL=hl.if_else(mt.col_idx == 0, 0.5, hl.missing(hl.tfloat64)))
    return mt


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        hl.init(tmp_dir=str(tmp_path / "hail_tmp"), quiet=True)

        mt = _build_mt()
        in_mt = tmp_path / "input.mt"
        mt.checkpoint(str(in_mt), overwrite=True)

        out_dir = tmp_path / "shards"
        args = argparse.Namespace(
            in_mt=str(in_mt),
            out_dir=str(out_dir),
            temp_dir=str(tmp_path / "spark_tmp"),
            shard_size=2,
            n_final_partitions=2,
            overwrite=True,
        )

        shard_main(args)

        shard_paths = sorted(out_dir.glob("shard_*.mt"))
        assert len(shard_paths) == 3

        counts = [hl.read_matrix_table(str(p)).count_cols() for p in shard_paths]
        assert counts == [2, 2, 1]

        print("smoke_test_shard_mt_by_samples: PASS")


if __name__ == "__main__":
    main()
