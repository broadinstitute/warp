#!/usr/bin/env python3
# v1.1.0

import argparse
import hail as hl


def build_default_metadata():
    """
    Metadata for the final exported VCF.

    AC/AF/AN/HC are recomputed.
    GQ and AVSAD are carried forward from the original INFO field.
    """
    return {
        "info": {
            "AC": {
                "Number": "A",
                "Type": "Integer",
                "Description": "Allele count in genotypes, for each ALT allele, in the same order as listed",
            },
            "AF": {
                "Number": "A",
                "Type": "Float",
                "Description": "Allele Frequency, for each ALT allele, in the same order as listed",
            },
            "AN": {
                "Number": "1",
                "Type": "Integer",
                "Description": "Total number of alleles in called genotypes",
            },
            "HC": {
                "Number": "R",
                "Type": "Integer",
                "Description": "Number of homozygotes per allele. One element per allele, including the reference.",
            },
            "AVSAD": {
                "Number": "1",
                "Type": "Float",
                "Description": "Mean sum of allelic depths. Proxies DP.",
            },
            "GQ": {
                "Number": "1",
                "Type": "Float",
                "Description": "Mean Genotype Quality",
            },
        },
        "filter": {
            "LowQual": {
                "Description": "Low quality score",
            },
            "ExcessHet": {
                "Description": "Excess heterozygotes",
            },
            "NO_HQ_GENOTYPES": {
                "Description": "No high-quality genotypes",
            },
        },
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Remove samples from a phased Hail MatrixTable, remove dead ALT alleles, "
            "recompute AC/AF/AN/HC, carry forward GQ/AVSAD, and export VCF."
        )
    )

    parser.add_argument(
        "--mt-path",
        required=True,
        help="Input Hail MatrixTable path.",
    )

    parser.add_argument(
        "--remove-samples-tsv",
        required=True,
        help="TSV containing sample/person IDs to remove.",
    )

    parser.add_argument(
        "--remove-id-col",
        default="research_id",
        help="Column in removal TSV matching mt.s. Default: research_id.",
    )

    parser.add_argument(
        "--out-vcf",
        required=True,
        help="Output VCF path, usually ending in .vcf.bgz.",
    )

    parser.add_argument(
        "--out-mt-path",
        default=None,
        help="Optional output MatrixTable path.",
    )

    parser.add_argument(
        "--metadata-vcf-or-header",
        default=None,
        help=(
            "Optional VCF or header file to use as metadata template. "
            "If omitted, a minimal metadata dictionary is used."
        ),
    )

    parser.add_argument(
        "--driver-memory",
        default="100g",
        help="Spark driver memory (legacy arg). Default: 100g.",
    )

    parser.add_argument(
        "--executor-memory",
        default="10g",
        help="Spark executor memory (legacy arg). Default: 10g.",
    )

    parser.add_argument(
        "--tmp-dir",
        default=None,
        help=(
            "Optional Hail temp directory (usually a GCS path for Dataproc runs), "
            "e.g. gs://bucket/path/tmp."
        ),
    )

    parser.add_argument(
        "--spark-local-dir",
        default=None,
        help=(
            "Optional Spark local scratch dir/path (often set to same GCS temp prefix in Dataproc wrapper)."
        ),
    )

    parser.add_argument(
        "--spark-driver-memory",
        default=None,
        help="Spark driver memory override, e.g. 60g.",
    )

    parser.add_argument(
        "--spark-driver-cores",
        type=int,
        default=None,
        help="Spark driver cores override, e.g. 8.",
    )

    parser.add_argument(
        "--spark-executor-memory",
        default=None,
        help="Spark executor memory override, e.g. 8g.",
    )

    parser.add_argument(
        "--spark-executor-cores",
        type=int,
        default=None,
        help="Spark executor cores override, e.g. 4.",
    )

    parser.add_argument(
        "--spark-task-max-failures",
        type=int,
        default=20,
        help="Spark task max failures. Default: 20.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing outputs.",
    )

    return parser.parse_args()


def build_spark_conf(args):
    spark_conf = {}

    driver_memory = args.spark_driver_memory or args.driver_memory
    executor_memory = args.spark_executor_memory or args.executor_memory

    if driver_memory:
        spark_conf["spark.driver.memory"] = driver_memory
    if executor_memory:
        spark_conf["spark.executor.memory"] = executor_memory
    if args.spark_driver_cores is not None:
        spark_conf["spark.driver.cores"] = str(args.spark_driver_cores)
    if args.spark_executor_cores is not None:
        spark_conf["spark.executor.cores"] = str(args.spark_executor_cores)
    if args.spark_local_dir:
        spark_conf["spark.local.dir"] = args.spark_local_dir
    if args.spark_task_max_failures is not None:
        spark_conf["spark.task.maxFailures"] = str(args.spark_task_max_failures)

    return spark_conf


def main():
    args = parse_args()

    spark_conf = build_spark_conf(args)

    print(f"[remove_phased_samples] mt_path={args.mt_path}")
    print(f"[remove_phased_samples] remove_samples_tsv={args.remove_samples_tsv}")
    print(f"[remove_phased_samples] out_vcf={args.out_vcf}")
    print(f"[remove_phased_samples] tmp_dir={args.tmp_dir}")
    print(f"[remove_phased_samples] spark_conf={spark_conf}")

    init_kwargs = {
        "default_reference": "GRCh38",
        "idempotent": True,
        "spark_conf": spark_conf,
    }
    if args.tmp_dir:
        init_kwargs["tmp_dir"] = args.tmp_dir

    hl.init(**init_kwargs)

    # ------------------------------------------------------------
    # 1. Read input MatrixTable
    # ------------------------------------------------------------
    mt = hl.read_matrix_table(args.mt_path)

    # ------------------------------------------------------------
    # 2. Import removal list
    # ------------------------------------------------------------
    # Force sample/person IDs to string so they match mt.s.
    remove_types = {
        args.remove_id_col: hl.tstr
    }

    person_ids_to_remove = hl.import_table(
        args.remove_samples_tsv,
        key=args.remove_id_col,
        types=remove_types,
        impute=False,
    )

    # ------------------------------------------------------------
    # 3. Remove samples
    # ------------------------------------------------------------
    mt = mt.filter_cols(
        ~hl.is_defined(person_ids_to_remove[mt.s])
    )

    # ------------------------------------------------------------
    # 4. Recompute call stats after sample removal
    # ------------------------------------------------------------
    # hl.agg.call_stats returns arrays that include REF at index 0.
    mt = mt.annotate_rows(
        cs=hl.agg.call_stats(mt.GT, mt.alleles)
    )

    # ------------------------------------------------------------
    # 5. Remove dead ALT alleles
    # ------------------------------------------------------------
    # Keep REF allele i == 0.
    # Keep ALT alleles only if recomputed AC > 0.
    mt_filtered = hl.filter_alleles(
        mt,
        lambda allele, i: (i == 0) | (mt.cs.AC[i] > 0),
    )

    # ------------------------------------------------------------
    # 6. Reindex GT after allele filtering
    # ------------------------------------------------------------
    # filter_alleles creates old_to_new and new_to_old mappings.
    # Preserve phasing.
    mt_filtered = mt_filtered.annotate_entries(
        GT=hl.if_else(
            hl.is_defined(mt_filtered.GT),
            hl.if_else(
                mt_filtered.GT.ploidy == 1,
                hl.call(
                    mt_filtered.old_to_new[mt_filtered.GT[0]],
                    phased=mt_filtered.GT.phased,
                ),
                hl.call(
                    mt_filtered.old_to_new[mt_filtered.GT[0]],
                    mt_filtered.old_to_new[mt_filtered.GT[1]],
                    phased=mt_filtered.GT.phased,
                ),
            ),
            hl.missing(hl.tcall),
        )
    )

    # ------------------------------------------------------------
    # 7. Remove rows that became REF-only
    # ------------------------------------------------------------
    mt_filtered = mt_filtered.filter_rows(
        hl.len(mt_filtered.alleles) > 1
    )

    # ------------------------------------------------------------
    # 8. Recompute final variant QC after allele removal / GT reindexing
    # ------------------------------------------------------------
    mt_filtered = hl.variant_qc(mt_filtered)

    # ------------------------------------------------------------
    # 9. Write final VCF-compatible INFO field
    # ------------------------------------------------------------
    # AC and AF are sliced [1:] because Hail includes REF at index 0,
    # while VCF INFO/AC and INFO/AF are ALT-indexed.
    #
    # HC is kept as Number=R, meaning one value per allele including REF.
    #
    # GQ and AVSAD are carried forward from the original INFO field.
    mt_filtered = mt_filtered.annotate_rows(
        info=mt_filtered.info.annotate(
            AC=mt_filtered.variant_qc.AC[1:],
            AF=mt_filtered.variant_qc.AF[1:],
            AN=mt_filtered.variant_qc.AN,
            HC=mt_filtered.variant_qc.homozygote_count,
            AVSAD=mt_filtered.info.AVSAD,
            GQ=mt_filtered.info.GQ,
        )
    )

    # ------------------------------------------------------------
    # 10. Drop temporary Hail fields
    # ------------------------------------------------------------
    fields_to_drop = [
        "cs",
        "variant_qc",
        "old_to_new",
        "new_to_old",
        "old_locus",
        "old_alleles",
    ]

    existing_fields_to_drop = [
        field for field in fields_to_drop
        if field in mt_filtered.row
    ]

    mt_filtered = mt_filtered.drop(*existing_fields_to_drop)

    # ------------------------------------------------------------
    # 11. Optionally write final MatrixTable
    # ------------------------------------------------------------
    if args.out_mt_path:
        mt_filtered.write(
            args.out_mt_path,
            overwrite=args.overwrite,
        )

    # ------------------------------------------------------------
    # 12. Export final VCF
    # ------------------------------------------------------------
    if args.metadata_vcf_or_header:
        metadata = hl.get_vcf_metadata(args.metadata_vcf_or_header)
    else:
        metadata = build_default_metadata()

    hl.export_vcf(
        mt_filtered,
        args.out_vcf,
        metadata=metadata,
        tabix=True,
    )


if __name__ == "__main__":
    main()