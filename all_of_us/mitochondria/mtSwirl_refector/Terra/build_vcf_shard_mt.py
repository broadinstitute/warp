#!/usr/bin/env python3
"""Build a multi-sample Hail MatrixTable shard from a shard manifest TSV.

Inputs
------
A shard manifest TSV with columns:
  - s   (sample ID)
  - vcf (path, typically gs://...)

Outputs
-------
Writes a Hail MatrixTable directory at --out-mt.

Design goals
------------
* Keep the same schema/formatting as the existing combine Step 3 output so that
  downstream merge stages are stable.
* Only do VCF import + basic canonicalization here.
* Do NOT apply covdb hom-ref imputation here (that happens once on the final MT).

"""

from __future__ import annotations

import argparse
import logging
import os

import hail as hl

from generate_mtdna_call_mt.merging_utils import vcf_merging


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Import one VCF shard manifest into a Hail MT")
    p.add_argument("--shard-tsv", required=True, help="Shard manifest TSV (columns: s, vcf)")
    p.add_argument("--out-mt", required=True, help="Output MT path (directory ending in .mt recommended)")
    p.add_argument("--temp-dir", required=True, help="Temp directory for Hail/Spark scratch")
    p.add_argument("--chunk-size", type=int, default=100, help="Chunk size for hierarchical merge (default: 100)")
    p.add_argument(
        "--include-extra-v2-fields",
        action="store_true",
        help="Include v2.1 extra entry fields (AD/OriginalSelfRefAlleles/SwappedFieldIDs)",
    )
    p.add_argument(
        "--n-final-partitions",
        type=int,
        default=128,
        help="Target partitions for the shard MT (default: 128)",
    )
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def main(args: argparse.Namespace) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(asctime)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
    )
    logger = logging.getLogger("build_vcf_shard_mt")

    hl.init(tmp_dir=args.temp_dir)

    logger.info("Reading shard manifest: %s", args.shard_tsv)
    ht = hl.import_table(args.shard_tsv, types={"s": hl.tstr, "vcf": hl.tstr}, impute=False)
    ht = ht.filter(hl.is_defined(ht.s) & hl.is_defined(ht.vcf) & (hl.len(ht.vcf) > 0))

    pairs = ht.select("s", "vcf")
    # Collect to driver: shard sizes are intentionally small (e.g., 2500)
    rows = pairs.aggregate(hl.agg.collect(hl.struct(s=pairs.s, vcf=pairs.vcf)))
    vcf_paths = {r.s: r.vcf for r in rows}

    logger.info("Shard contains %d samples", len(vcf_paths))

    if len(vcf_paths) == 0:
        raise ValueError("Shard manifest produced 0 usable (s, vcf) rows")

    out = args.out_mt
    if hl.hadoop_exists(f"{out}/_SUCCESS") and not args.overwrite:
        logger.info("Output exists, reading: %s", out)
        _ = hl.read_matrix_table(out)
        return

    logger.info("Importing and merging VCFs for shard...")
    combined_mt, _meta = vcf_merging(
        vcf_paths=vcf_paths,
        temp_dir=args.temp_dir,
        logger=logger,
        n_final_partitions=args.n_final_partitions,
        chunk_size=args.chunk_size,
        include_extra_v2_fields=args.include_extra_v2_fields,
        num_merges=1,
        single_sample=True,
    )

    logger.info("Checkpointing shard MT: %s", out)
    combined_mt = combined_mt.naive_coalesce(args.n_final_partitions).checkpoint(out, overwrite=args.overwrite)


if __name__ == "__main__":
    main(_parse_args())
