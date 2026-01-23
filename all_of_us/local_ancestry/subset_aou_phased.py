#!/usr/bin/env python3
"""
Subset AoU per-chrom phased VCFs to a provided sample list, optionally checkpoint,
export a bgzipped VCF.

Designed to be called from WDL:
- all input/output paths are CLI args
- supports GCS paths (gs://...) for inputs/outputs

Notes
- hl.export_vcf writes block-gzipped VCF when the output ends with .bgz/.gz.
- If you pass --checkpoint_path, the subset MT will be checkpointed before export.
"""

import argparse
import logging
from typing import Optional, List

import hail as hl


# -----------------------------
# Hail init
# -----------------------------
def hail_init(
    executor_memory: str,
    executor_cores: str,
    driver_cores: str,
    driver_memory: str,
    reference_genome: str,
    *,
    min_partitions: Optional[int] = None,
    tmp_dir: Optional[str] = None,
    quiet: bool = False,
) -> None:
    spark_conf = {
        "spark.executor.memory": executor_memory,
        "spark.executor.cores": executor_cores,
        "spark.driver.memory": driver_memory,
        "spark.driver.cores": driver_cores,
    }
    if min_partitions is not None:
        spark_conf["spark.sql.shuffle.partitions"] = str(min_partitions)

    hl.init(
        default_reference=reference_genome,
        idempotent=True,
        spark_conf=spark_conf,
        tmp_dir=tmp_dir,
        quiet=quiet,
        skip_logging_configuration=False,
    )


# -----------------------------
# VCF subsetting helpers
# -----------------------------
def vcf_path_for_chrom(phased_vcf_gcs_dir: str, chrom: str) -> str:
    phased_vcf_gcs_dir = phased_vcf_gcs_dir.rstrip("/")
    return f"{phased_vcf_gcs_dir}/{chrom}.aou.v9.phased.vcf.gz"


def import_phased_vcf_chr(
    phased_vcf_gcs_dir: str,
    chrom: str,
    *,
    reference_genome: str = "GRCh38",
    force_bgz: bool = True,
    min_partitions: Optional[int] = None,
    entry_fields: Optional[List[str]] = None,
) -> hl.MatrixTable:
    path = vcf_path_for_chrom(phased_vcf_gcs_dir, chrom)

    mt = hl.import_vcf(
        path,
        reference_genome=reference_genome,
        force_bgz=force_bgz,
        min_partitions=min_partitions,
    )

    if entry_fields:
        mt = mt.select_entries(**{f: mt[f] for f in entry_fields})

    return mt


def load_samples_table(
    samples_tsv: str,
    sample_id_column: str,
    *,
    delimiter: str = "\t",
    comment: Optional[str] = None,
    no_header: bool = False,
) -> hl.Table:
    kwargs = dict(
        impute=True,
        delimiter=delimiter,
        no_header=no_header,
    )
    if comment is not None:
        kwargs["comment"] = comment

    ht = hl.import_table(samples_tsv, **kwargs)

    if sample_id_column not in ht.row:
        raise ValueError(
            f"Column '{sample_id_column}' not found in {samples_tsv}. "
            f"Found columns: {list(ht.row)}"
        )

    return ht.annotate(_sid=hl.str(ht[sample_id_column])).key_by("_sid")


def subset_mt_to_samples(mt: hl.MatrixTable, samples_ht: hl.Table, *, keep: bool = True) -> hl.MatrixTable:
    return mt.semi_join_cols(samples_ht) if keep else mt.anti_join_cols(samples_ht)


def subset_phased_vcf_chr_to_samples(
    phased_vcf_gcs_dir: str,
    chrom: str,
    samples_tsv: str,
    sample_id_column: str,
    *,
    reference_genome: str = "GRCh38",
    force_bgz: bool = True,
    min_partitions: Optional[int] = None,
    entry_fields: Optional[List[str]] = None,
) -> hl.MatrixTable:
    mt = import_phased_vcf_chr(
        phased_vcf_gcs_dir=phased_vcf_gcs_dir,
        chrom=chrom,
        reference_genome=reference_genome,
        force_bgz=force_bgz,
        min_partitions=min_partitions,
        entry_fields=entry_fields,
    )
    samples_ht = load_samples_table(samples_tsv, sample_id_column)
    return subset_mt_to_samples(mt, samples_ht, keep=True)


# -----------------------------
# CLI
# -----------------------------
def parse_arguments() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Subset per-chrom phased VCF to samples; optional checkpoint; export VCF."
    )

    # Hail/Spark
    p.add_argument("--executor_memory", required=True, help="e.g. 8g")
    p.add_argument("--executor_cores", required=True, help="e.g. 4")
    p.add_argument("--driver_memory", required=True, help="e.g. 120g")
    p.add_argument("--driver_cores", required=True, help="e.g. 8")
    p.add_argument("--reference_genome", default="GRCh38", help='e.g. "GRCh38"')
    p.add_argument("--min_partitions", type=int, default=None, help="min_partitions for import_vcf / tuning")
    p.add_argument("--tmp_dir", default=None, help="Hail tmp_dir (recommend a gs://... path in Terra)")
    p.add_argument("--quiet", action="store_true", help="Reduce Hail logging verbosity")

    # VCF subsetting inputs
    p.add_argument("--phased_vcf_gcs_dir", required=True, help='e.g. "gs://prod-drc-broad/aou_phasing/v9"')
    p.add_argument("--chrom", required=True, help='e.g. "chr21"')
    p.add_argument("--samples_tsv", required=True, help="TSV containing sample IDs (gs://... or local)")
    p.add_argument("--sample_id_column", required=True, help='Column name in TSV, e.g. "research_id"')

    # Optional VCF import/export tuning
    p.add_argument(
        "--entry_fields",
        default="GT",
        help='Comma-separated entry fields to keep (default: "GT"). Use "" to keep all.',
    )
    p.add_argument(
        "--force_bgz",
        action="store_true",
        help="Pass force_bgz=True to hl.import_vcf (useful if extension detection is unreliable).",
    )
    p.add_argument("--output_vcf", required=True, help="Where to write the subset VCF (gs://.../*.vcf.bgz)")

    # Optional checkpoint
    p.add_argument("--checkpoint_path", default=None, help="Optional: checkpoint MT to this path (gs://.../*.mt)")
    p.add_argument("--checkpoint_overwrite", action="store_true", help="Overwrite checkpoint if it exists")

    # Generic
    p.add_argument("--log_level", default="INFO", help="DEBUG/INFO/WARNING/ERROR")
    return p.parse_args()


def main() -> None:
    args = parse_arguments()
    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO))
    logger = logging.getLogger("subset_vcf_hail")

    hail_init(
        args.executor_memory,
        args.executor_cores,
        args.driver_cores,
        args.driver_memory,
        args.reference_genome,
        min_partitions=args.min_partitions,
        tmp_dir=args.tmp_dir,
        quiet=args.quiet,
    )

    # Parse entry_fields
    entry_fields = None
    if args.entry_fields is not None:
        ef = [x.strip() for x in args.entry_fields.split(",") if x.strip() != ""]
        entry_fields = ef if len(ef) > 0 else None

    logger.info(
        "Subsetting VCF: dir=%s chrom=%s samples=%s col=%s",
        args.phased_vcf_gcs_dir,
        args.chrom,
        args.samples_tsv,
        args.sample_id_column,
    )

    mt_sub = subset_phased_vcf_chr_to_samples(
        phased_vcf_gcs_dir=args.phased_vcf_gcs_dir,
        chrom=args.chrom,
        samples_tsv=args.samples_tsv,
        sample_id_column=args.sample_id_column,
        reference_genome=args.reference_genome,
        force_bgz=args.force_bgz,
        min_partitions=args.min_partitions,
        entry_fields=entry_fields,
    )

    # Optional checkpoint
    mt_to_export = mt_sub
    if args.checkpoint_path:
        logger.info("Checkpointing MT to %s", args.checkpoint_path)
        mt_to_export = mt_sub.checkpoint(args.checkpoint_path, overwrite=args.checkpoint_overwrite)

    logger.info("Exporting VCF to %s", args.output_vcf)
    hl.export_vcf(mt_to_export, args.output_vcf)

    logger.info("Done.")


if __name__ == "__main__":
    main()
