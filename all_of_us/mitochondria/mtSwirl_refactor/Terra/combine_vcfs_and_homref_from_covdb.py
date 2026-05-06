#!/usr/bin/env python3
"""Combine per-sample mitochondrial VCFs into a single Hail MatrixTable and
impute homoplasmic reference calls using the v2 HDF5 coverage DB.

This is the v2 replacement for `Terra/combine_vcfs.py` + `merging_utils.determine_hom_refs`.

Key contracts (must match v1 output semantics):
- Input VCF paths come from a TSV column (often gs:// URIs) and will not be localized by WDL.
- Hom-ref/DP imputation matches `determine_hom_refs` exactly:
    * If HL is missing, set DP from coverage.
    * If HL is missing and DP > minimum_homref_coverage => HL=0.0, FT={"PASS"}.
    * If HL is missing and DP <= minimum_homref_coverage => DP is set back to missing.

Primary output:
- A Hail MatrixTable written to `output_bucket/<file_name>.mt`.

Notes on scaling:
- This keeps the v1 hierarchical MT merge logic.
- It avoids building a coverage MatrixTable by sourcing coverage from `coverage.h5`.

"""

from __future__ import annotations

import argparse
import logging
import os

import hail as hl

from generate_mtdna_call_mt.merging_utils import (
    collect_vcf_paths,
    determine_hom_refs_from_covdb,
    vcf_merging_and_processing,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("combine_vcfs_and_homref_from_covdb")


def main(args: argparse.Namespace) -> None:
    participant_data = args.input_tsv
    coverage_h5_path = args.coverage_h5_path
    output_bucket = args.output_bucket
    temp_dir = args.temp_dir
    participants_to_subset = None if args.participants_to_subset is None else f"gs://{args.participants_to_subset}"
    chunk_size = args.chunk_size
    artifact_prone_sites_path = args.artifact_prone_sites_path
    artifact_prone_sites_reference = args.artifact_prone_sites_reference
    include_extra_v2_fields = args.include_extra_v2_fields
    file_name = args.file_name
    minimum_homref_coverage = args.minimum_homref_coverage
    num_merges = args.split_merging
    vcf_col_name = args.vcf_col_name

    hl.init(tmp_dir=temp_dir)

    if int(hl.version().split("-")[0].split(".")[2]) >= 75:
        logger.info("Setting hail flag to avoid array index out of bounds error...")
        hl._set_flags(no_whole_stage_codegen="1")

    logger.info("Collecting VCF paths for samples to subset...")
    vcf_paths = collect_vcf_paths(
        participant_data, vcf_col_name, participants_to_subset, single_sample=True
    )
    logger.info("Found %d VCF paths", len(vcf_paths))

    if args.append_to_existing is not None:
        logger.info("Will search for existing mt to join new data into.")

    combined_mt, meta = vcf_merging_and_processing(
        vcf_paths=vcf_paths,
        coverage_mt_path=None,  # v2: coverage MT not used
        include_extra_v2_fields=include_extra_v2_fields,
        single_sample=True,
        old_mt_path=args.append_to_existing,
        artifact_prone_sites_path=artifact_prone_sites_path,
        artifact_prone_sites_reference=artifact_prone_sites_reference,
        minimum_homref_coverage=minimum_homref_coverage,
        logger=logger,
        chunk_size=chunk_size,
        num_merges=num_merges,
        n_final_partitions=args.n_final_partitions,
        output_bucket=output_bucket,
        temp_dir=temp_dir,
        overwrite=args.overwrite,
        # v2-only:
        coverage_h5_path=coverage_h5_path,
        homref_position_block_size=args.homref_position_block_size,
    )

    logger.info("Writing combined MT...")
    out_mt = os.path.join(output_bucket, f"{file_name}.mt")
    combined_mt = combined_mt.repartition(args.n_final_partitions).checkpoint(out_mt, overwrite=args.overwrite)

    # Optional: write a trimmed entries TSV for debugging/side products (kept from v1).
    if args.write_entries_tsv:
        logger.info("Writing trimmed variants table...")
        out_tsv = os.path.join(output_bucket, f"{file_name}.tsv.bgz")
        ht_for_tsv = combined_mt.entries()
        ht_for_tsv = ht_for_tsv.filter(hl.is_missing(ht_for_tsv.HL) | (ht_for_tsv.HL > 0))
        ht_for_tsv.repartition(50).export(out_tsv)

    logger.info("Done. MT written to %s", out_mt)


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description=(
            "Combine individual mitochondria VCF files into one MatrixTable, "
            "determine homoplasmic reference sites using coverage.h5, and apply artifact_prone_site filter"
        )
    )
    p.add_argument(
        "-i",
        "--input-tsv",
        help="Input tsv with paths to vcf file.",
        required=True,
    )
    p.add_argument(
        "--coverage-h5-path",
        help="Path to coverage.h5 produced by build_coverage_db (local path after extraction/localization).",
        required=True,
    )
    p.add_argument(
        "-v",
        "--vcf-col-name",
        help="Name of column in participant data file that contains the path to the VCF output by Mutect2",
        required=True,
    )
    p.add_argument(
        "-a",
        "--artifact-prone-sites-path",
        help="Path to BED file of artifact-prone sites to flag in the FILTER column",
        required=True,
    )
    p.add_argument(
        "-a-ref",
        "--artifact-prone-sites-reference",
        help="Reference genome for artifact-prone sites BED file. If specified can be GRCh37 or GRCh38",
        default="default",
        type=str,
    )
    p.add_argument(
        "-o",
        "--output-bucket",
        help="Path to bucket to which results should be written",
        required=True,
    )
    p.add_argument(
        "-t",
        "--temp-dir",
        help="Temporary directory to use for intermediate outputs",
        required=True,
    )
    p.add_argument(
        "-f",
        "--file-name",
        help="File name to use for output MT (will be used for the .mt output)",
        required=True,
    )
    p.add_argument(
        "-s",
        "--participants-to-subset",
        help="Path to txt file of participant_ids to which the data should be subset (file should contain header (named 'participant') and one line for each participant_id matching the 'entity:participant_id's supplied in Terra",
    )
    p.add_argument(
        "--minimum-homref-coverage",
        help="Minimum depth of coverage required to call a genotype homoplasmic reference rather than missing",
        type=int,
        default=100,
    )
    p.add_argument(
        "--chunk-size",
        help="Chunk size to use for combining VCFs (the number of individual VCFs that should be combined at a time)",
        type=int,
        default=100,
    )
    p.add_argument(
        "--homref-position-block-size",
        help="Number of mitochondrial positions to process per block when applying hom-ref/DP imputation from coverage.h5",
        type=int,
        default=1024,
    )
    p.add_argument("--overwrite", help="Overwrites existing files", action="store_true")
    p.add_argument(
        "--include-extra-v2-fields",
        help=(
            "Loads in extra fields from MitochondriaPipeline v2.1, specifically AD, "
            "OriginalSelfRefAlleles, and SwappedFieldIDs. If missing will fill with missing."
        ),
        action="store_true",
    )
    p.add_argument(
        "--n-final-partitions", type=int, default=1000, help="Number of partitions for final mt."
    )
    p.add_argument(
        "--split-merging",
        type=int,
        default=1,
        help=(
            "Will split the merging into this many jobs which will be merged at the end. "
            "Uses the same order each time such that if it fails we can read from previous files."
        ),
    )
    p.add_argument(
        "--append-to-existing",
        type=str,
        help="Optional: specify absolute path to existing MatrixTable to merge this dataset into.",
    )
    p.add_argument(
        "--write-entries-tsv",
        action="store_true",
        help="Optional: export a trimmed entries TSV (kept from v1 behavior).",
    )

    args = p.parse_args()
    main(args)
