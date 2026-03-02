#!/usr/bin/env python3
"""Finalize a merged callset MT by applying covdb-based hom-ref/DP and artifact-prone filter.

This step runs once, after all shard merges, to preserve exact cohort-wide behavior.

Inputs
------
* --in-mt: merged MT (from reduce-tree merging)
* --coverage-h5-path: local path to coverage.h5 (localized/untarred by WDL)

Output
------
* --out-mt: finalized MT checkpoint

"""

from __future__ import annotations

import argparse
import logging

import hail as hl

from generate_mtdna_call_mt.merging_utils import (
    apply_mito_artifact_filter,
    determine_hom_refs_from_covdb,
    remove_genotype_filters,
)


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Apply covdb homref/DP + artifact filter to a merged MT")
    p.add_argument("--in-mt", required=True, help="Input merged MT")
    p.add_argument("--coverage-h5-path", required=True, help="Local path to coverage.h5")
    p.add_argument("--out-mt", required=True, help="Output finalized MT")
    p.add_argument("--temp-dir", required=True, help="Temp directory for Hail/Spark scratch")
    p.add_argument("--minimum-homref-coverage", type=int, default=100)
    p.add_argument("--homref-position-block-size", type=int, default=1024)
    p.add_argument("--covdb-checkpoint-interval-blocks", type=int, default=1)
    p.add_argument(
        "--artifact-prone-sites-path",
        default="gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.bed",
    )
    p.add_argument("--artifact-prone-sites-reference", default="default")
    p.add_argument("--n-final-partitions", type=int, default=1000)
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def main(args: argparse.Namespace) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(asctime)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
    )
    logger = logging.getLogger("finalize_mt_with_covdb")

    hl.init(tmp_dir=args.temp_dir)

    if hl.hadoop_exists(f"{args.out_mt}/_SUCCESS") and not args.overwrite:
        logger.info("Output exists, reading: %s", args.out_mt)
        _ = hl.read_matrix_table(args.out_mt)
        return

    logger.info("Reading merged MT: %s", args.in_mt)
    mt = hl.read_matrix_table(args.in_mt)

    logger.info("Removing select sample-level filters...")
    mt = remove_genotype_filters(mt)

    logger.info("Applying covdb-based homref/DP imputation...")
    mt = determine_hom_refs_from_covdb(
        mt,
        coverage_h5_path=args.coverage_h5_path,
        minimum_homref_coverage=args.minimum_homref_coverage,
        position_block_size=args.homref_position_block_size,
        checkpoint_interval_blocks=args.covdb_checkpoint_interval_blocks,
        logger=logger,
    )

    logger.info("Applying artifact_prone_site filter...")
    mt = apply_mito_artifact_filter(
        mt, args.artifact_prone_sites_path, args.artifact_prone_sites_reference
    )

    # Defensive guard: ensure row keys are fully defined before checkpoint/write.
    # In Hail 0.2.127, missing keys can trigger JVM-side NPEs during RVD key coercion.
    logger.info("Filtering rows with missing keys before checkpoint...")
    mt = mt.filter_rows(hl.is_defined(mt.locus) & hl.is_defined(mt.alleles))

    logger.info("Checkpointing finalized MT: %s", args.out_mt)
    mt = mt.repartition(args.n_final_partitions)
    mt.checkpoint(args.out_mt, overwrite=args.overwrite, stage_locally=True)


if __name__ == "__main__":
    main(_parse_args())
