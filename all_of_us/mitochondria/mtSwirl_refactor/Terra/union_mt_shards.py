#!/usr/bin/env python3
"""Union MT shards by columns with row-key alignment.

Inputs
------
* --mt-list-tsv: TSV with header 'mt_path' and one MT path per line.

Output
------
* --out-mt: merged MT directory.
"""

from __future__ import annotations

import argparse
import logging

import hail as hl


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Union MT shards by columns")
    p.add_argument("--mt-list-tsv", required=True, help="TSV with header mt_path")
    p.add_argument("--out-mt", required=True, help="Output MT path")
    p.add_argument("--temp-dir", required=True, help="Temp directory for Hail/Spark scratch")
    p.add_argument(
        "--n-final-partitions",
        type=int,
        default=1000,
        help="Target partitions for the merged MT",
    )
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def main(args: argparse.Namespace) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(asctime)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
    )
    logger = logging.getLogger("union_mt_shards")

    hl.init(tmp_dir=args.temp_dir)

    if hl.hadoop_exists(f"{args.out_mt}/_SUCCESS") and not args.overwrite:
        logger.info("Output exists, reading: %s", args.out_mt)
        _ = hl.read_matrix_table(args.out_mt)
        return

    ht = hl.import_table(args.mt_list_tsv, impute=False, types={"mt_path": hl.tstr})
    paths = ht.aggregate(hl.agg.collect(ht.mt_path))
    paths = [p for p in paths if p is not None and len(p) > 0]

    if len(paths) == 0:
        raise ValueError("mt-list-tsv contained 0 mt_path entries")

    logger.info("Unioning %d MT shards", len(paths))
    mts = [hl.read_matrix_table(p) for p in paths]

    merged = mts[0]
    for mt in mts[1:]:
        merged = merged.union_cols(mt, row_join_type="outer")

    logger.info("Checkpoint merged MT: %s", args.out_mt)
    merged = merged.repartition(args.n_final_partitions)
    merged.checkpoint(args.out_mt, overwrite=args.overwrite)


if __name__ == "__main__":
    main(_parse_args())
