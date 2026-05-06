#!/usr/bin/env python3
"""Merge multiple Hail MatrixTable shards into a single MatrixTable.

This is intended for a WDL reduce-tree pattern (fan-in merging).

Inputs
------
* --mt-list-tsv: TSV with a header 'mt_path' and one MT path per line.

Output
------
* Writes the merged MatrixTable to --out-mt.

Notes
-----
We reuse the existing v1 multi-way zip merge implementation (`multi_way_union_mts`)
so we keep stable formatting.

"""

from __future__ import annotations

import argparse
import logging

import hail as hl

from generate_mtdna_call_mt.merging_utils import multi_way_union_mts


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Merge a list of MTs into one MT")
    p.add_argument("--mt-list-tsv", required=True, help="TSV with header mt_path")
    p.add_argument("--out-mt", required=True, help="Output MT path")
    p.add_argument("--temp-dir", required=True, help="Temp directory for Hail/Spark scratch")
    p.add_argument("--chunk-size", type=int, default=10, help="Fan-in chunk size per merge stage (default: 10)")
    p.add_argument(
        "--min-partitions",
        type=int,
        default=256,
        help="Minimum partitions to use during merge checkpoints (default: 256)",
    )
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def main(args: argparse.Namespace) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(asctime)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
    )
    logger = logging.getLogger("merge_mt_shards")

    hl.init(tmp_dir=args.temp_dir)

    ht = hl.import_table(args.mt_list_tsv, impute=False, types={"mt_path": hl.tstr})
    paths = ht.aggregate(hl.agg.collect(ht.mt_path))
    paths = [p for p in paths if p is not None and len(p) > 0]

    if len(paths) == 0:
        raise ValueError("mt-list-tsv contained 0 mt_path entries")

    if hl.hadoop_exists(f"{args.out_mt}/_SUCCESS") and not args.overwrite:
        logger.info("Output exists, reading: %s", args.out_mt)
        _ = hl.read_matrix_table(args.out_mt)
        return

    logger.info("Merging %d MTs", len(paths))
    mts = [hl.read_matrix_table(p) for p in paths]

    merged = multi_way_union_mts(
        mts,
        temp_dir=args.temp_dir,
        chunk_size=args.chunk_size,
        min_partitions=args.min_partitions,
        check_from_disk=False,
        prefix="merge_round/",
    )

    logger.info("Checkpoint merged MT: %s", args.out_mt)
    merged.checkpoint(args.out_mt, overwrite=args.overwrite)


if __name__ == "__main__":
    main(_parse_args())
