#!/usr/bin/env python3
"""Shard a MatrixTable by columns (samples).

This reads a single MT and writes multiple MT shards, each containing a subset
of samples. Shards are written to --out-dir as shard_XXXXX.mt directories.
"""

from __future__ import annotations

import argparse
import logging
import math
import os
from typing import List

import hail as hl


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Shard a MatrixTable by columns")
    p.add_argument("--in-mt", required=True, help="Input MT path")
    p.add_argument("--out-dir", required=True, help="Output directory for shard MTs")
    p.add_argument("--temp-dir", required=True, help="Temp directory for Hail/Spark scratch")
    p.add_argument("--shard-size", type=int, default=25000, help="Samples per shard")
    p.add_argument(
        "--n-final-partitions",
        type=int,
        default=256,
        help="Target partitions per shard",
    )
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def _chunk_list(items: List[str], chunk_size: int) -> List[List[str]]:
    return [items[i : i + chunk_size] for i in range(0, len(items), chunk_size)]


def main(args: argparse.Namespace) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(asctime)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
    )
    logger = logging.getLogger("shard_mt_by_samples")

    if args.shard_size <= 0:
        raise ValueError("shard-size must be > 0")

    os.makedirs(args.out_dir, exist_ok=True)

    hl.init(tmp_dir=args.temp_dir)

    logger.info("Reading MT: %s", args.in_mt)
    mt = hl.read_matrix_table(args.in_mt)

    mt = mt.add_col_index(name="__col_idx")
    col_rows = mt.cols().select(s=mt.s, idx=mt.__col_idx).collect()
    col_rows.sort(key=lambda r: int(r.idx))
    samples = [r.s for r in col_rows]

    if not samples:
        raise ValueError("MT has 0 columns; cannot shard")

    n_shards = int(math.ceil(len(samples) / args.shard_size))
    logger.info("Sharding %d samples into %d shards", len(samples), n_shards)

    for shard_id, shard_samples in enumerate(_chunk_list(samples, args.shard_size)):
        shard_name = f"shard_{shard_id:05d}.mt"
        out_mt = os.path.join(args.out_dir, shard_name)

        if hl.hadoop_exists(f"{out_mt}/_SUCCESS") and not args.overwrite:
            logger.info("Shard exists, skipping: %s", out_mt)
            continue

        sample_set = hl.literal(set(shard_samples))
        mt_shard = mt.filter_cols(sample_set.contains(mt.s)).drop("__col_idx")
        mt_shard = mt_shard.naive_coalesce(args.n_final_partitions)

        logger.info(
            "Writing shard %d/%d (%d samples) -> %s",
            shard_id + 1,
            n_shards,
            len(shard_samples),
            out_mt,
        )
        mt_shard.checkpoint(out_mt, overwrite=args.overwrite)


if __name__ == "__main__":
    main(_parse_args())
