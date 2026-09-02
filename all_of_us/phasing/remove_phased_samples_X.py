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
            "Remove samples from a phased chrX Hail MatrixTable, temporarily collapse "
            "male non-PAR genotypes to haploid for metric recalculation, remove dead "
            "ALT alleles, recompute AC/AF/AN/HC, restore homozygous diploid male GT "
            "for export, and carry forward GQ/AVSAD."
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
        "--sex-tsv",
        required=True,
        help="TSV containing sex info for every sample remaining in the MT.",
    )

    parser.add_argument(
        "--sex-id-col",
        default="research_id",
        help="Column in sex TSV matching mt.s. Default: research_id.",
    )

    parser.add_argument(
        "--sex-col",
        default="sex_at_birth",
        help="Column in sex TSV containing the sex value. Default: sex_at_birth.",
    )

    parser.add_argument(
        "--male-values",
        default="M,Male,MALE",
        help="Comma-separated list of sex values that indicate male.",
    )

    parser.add_argument(
        "--test-2kb-region",
        action="store_true",
        help=(
            "If set, filter rows to a small chrX interval before processing. "
            "Useful for quick test runs."
        ),
    )

    parser.add_argument(
        "--test-interval",
        default="chrX:1-2000",
        help=(
            "Interval used when --test-2kb-region is set. "
            "Default: chrX:1-2000."
        ),
    )

    parser.add_argument(
        "--test-male-count-tsv",
        default=None,
        help=(
            "Optional output TSV path for test mode male count summary. "
            "Writes one row with test_interval and males_marked."
        ),
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
        default="60g",
        help="Spark driver memory (legacy arg). Default: 60g.",
    )

    parser.add_argument(
        "--executor-memory",
        default="8g",
        help="Spark executor memory (legacy arg). Default: 8g.",
    )

    parser.add_argument(
        "--temp-bucket",
        default=None,
        help=(
            "GCS bucket/prefix used to derive AoU-style temp dir '<temp_bucket>/Stage_1/temp'. "
            "Used when --tmp-dir is not provided."
        ),
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
        default=32,
        help="Spark driver cores override. Default: 32.",
    )

    parser.add_argument(
        "--spark-executor-memory",
        default=None,
        help="Spark executor memory override, e.g. 8g.",
    )

    parser.add_argument(
        "--spark-executor-cores",
        type=int,
        default=4,
        help="Spark executor cores override. Default: 4.",
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
    effective_tmp_dir = get_effective_tmp_dir(args)
    spark_local_dir = args.spark_local_dir or effective_tmp_dir
    if spark_local_dir:
        spark_conf["spark.local.dir"] = spark_local_dir
    if args.spark_task_max_failures is not None:
        spark_conf["spark.task.maxFailures"] = str(args.spark_task_max_failures)

    return spark_conf


def get_effective_tmp_dir(args):
    if args.tmp_dir:
        return args.tmp_dir
    if args.temp_bucket:
        return f"{args.temp_bucket.rstrip('/')}/Stage_1/temp"
    return None


def main():
    args = parse_args()
    effective_tmp_dir = get_effective_tmp_dir(args)

    spark_conf = build_spark_conf(args)
    male_values = {value.strip().lower() for value in args.male_values.split(",") if value.strip()}

    print(f"[remove_phased_samples_X] mt_path={args.mt_path}")
    print(f"[remove_phased_samples_X] remove_samples_tsv={args.remove_samples_tsv}")
    print(f"[remove_phased_samples_X] sex_tsv={args.sex_tsv}")
    print(f"[remove_phased_samples_X] test_2kb_region={args.test_2kb_region}")
    if args.test_2kb_region:
        print(f"[remove_phased_samples_X] test_interval={args.test_interval}")
        if args.test_male_count_tsv:
            print(f"[remove_phased_samples_X] test_male_count_tsv={args.test_male_count_tsv}")
    print(f"[remove_phased_samples_X] out_vcf={args.out_vcf}")
    print(f"[remove_phased_samples_X] tmp_dir={effective_tmp_dir}")
    print(f"[remove_phased_samples_X] spark_conf={spark_conf}")

    init_kwargs = {
        "default_reference": "GRCh38",
        "idempotent": True,
        "spark_conf": spark_conf,
    }
    if effective_tmp_dir:
        init_kwargs["tmp_dir"] = effective_tmp_dir

    hl.init(**init_kwargs)

    # ------------------------------------------------------------
    # 1. Read input MatrixTable
    # ------------------------------------------------------------
    mt = hl.read_matrix_table(args.mt_path)

    # ------------------------------------------------------------
    # 1b. Optional test filter to a small chrX interval
    # ------------------------------------------------------------
    if args.test_2kb_region:
        test_interval = hl.parse_locus_interval(
            args.test_interval,
            reference_genome="GRCh38",
        )
        mt = hl.filter_intervals(mt, [test_interval])

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
    # 4. Annotate sex and convert male non-PAR chrX calls to haploid
    # ------------------------------------------------------------
    sex_types = {
        args.sex_id_col: hl.tstr,
    }

    sex_ht = hl.import_table(
        args.sex_tsv,
        key=args.sex_id_col,
        types=sex_types,
        impute=False,
    )

    sex_ht = sex_ht.annotate(
        is_male=hl.if_else(
            hl.is_defined(sex_ht[args.sex_col]),
            hl.literal(male_values).contains(hl.str(sex_ht[args.sex_col]).lower()),
            hl.missing(hl.tbool),
        )
    )

    mt = mt.annotate_cols(
        is_male=sex_ht[mt.s].is_male
    )

    n_missing_sex = mt.aggregate_cols(
        hl.agg.count_where(~hl.is_defined(mt.is_male))
    )
    if n_missing_sex > 0:
        raise ValueError(
            f"{n_missing_sex} sample(s) in the MT (after removal) have no matching "
            f"entry in --sex-tsv."
        )

    if args.test_2kb_region:
        n_males = mt.aggregate_cols(
            hl.agg.count_where(mt.is_male)
        )
        print(f"[remove_phased_samples_X] test_mode_males_marked={n_males}")

        if args.test_male_count_tsv:
            male_count_ht = hl.utils.range_table(1, n_partitions=1).key_by().select(
                test_interval=args.test_interval,
                males_marked=n_males,
            )
            male_count_ht.export(args.test_male_count_tsv, header=True)
            print(f"[remove_phased_samples_X] wrote_test_male_count_tsv={args.test_male_count_tsv}")

    mt = mt.annotate_rows(
        in_x_nonpar=mt.locus.in_x_nonpar()
    )

    mt = mt.annotate_entries(
        GT=hl.case()
        .when(~(mt.is_male & mt.in_x_nonpar), mt.GT)
        .when(~hl.is_defined(mt.GT), mt.GT)
        .when(mt.GT.is_het(), hl.missing(hl.tcall))
        .default(hl.call(mt.GT[0], phased=mt.GT.phased))
    )

    # ------------------------------------------------------------
    # 5. Recompute call stats after sample removal
    # ------------------------------------------------------------
    # hl.agg.call_stats returns arrays that include REF at index 0.
    mt = mt.annotate_rows(
        cs=hl.agg.call_stats(mt.GT, mt.alleles)
    )

    # ------------------------------------------------------------
    # 6. Remove dead ALT alleles
    # ------------------------------------------------------------
    # Keep REF allele i == 0.
    # Keep ALT alleles only if recomputed AC > 0.
    mt_filtered = hl.filter_alleles(
        mt,
        lambda allele, i: (i == 0) | (mt.cs.AC[i] > 0),
    )

    # ------------------------------------------------------------
    # 7. Reindex GT after allele filtering
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
    # 8. Remove rows that became REF-only
    # ------------------------------------------------------------
    mt_filtered = mt_filtered.filter_rows(
        hl.len(mt_filtered.alleles) > 1
    )

    # ------------------------------------------------------------
    # 9. Recompute final variant QC after allele removal / GT reindexing
    # ------------------------------------------------------------
    mt_filtered = hl.variant_qc(mt_filtered)

    # ------------------------------------------------------------
    # 10. Write final VCF-compatible INFO field
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
    # 11. Restore homozygous diploid GT for male non-PAR samples
    # ------------------------------------------------------------
    mt_filtered = mt_filtered.annotate_entries(
        GT=hl.if_else(
            mt_filtered.is_male & mt_filtered.in_x_nonpar,
            hl.if_else(
                hl.is_defined(mt_filtered.GT),
                hl.if_else(
                    mt_filtered.GT.ploidy == 1,
                    hl.call(
                        mt_filtered.GT[0],
                        mt_filtered.GT[0],
                        phased=mt_filtered.GT.phased,
                    ),
                    mt_filtered.GT,
                ),
                hl.missing(hl.tcall),
            ),
            mt_filtered.GT,
        )
    )

    # ------------------------------------------------------------
    # 12. Drop temporary Hail fields
    # ------------------------------------------------------------
    fields_to_drop = [
        "cs",
        "variant_qc",
        "is_male",
        "in_x_nonpar",
        "old_to_new",
        "new_to_old",
        "old_locus",
        "old_alleles",
    ]

    existing_row_fields_to_drop = [
        field for field in fields_to_drop
        if field in mt_filtered.row
    ]

    existing_col_fields_to_drop = [
        field for field in fields_to_drop
        if field in mt_filtered.col and field not in existing_row_fields_to_drop
    ]

    mt_filtered = mt_filtered.drop(*existing_row_fields_to_drop, *existing_col_fields_to_drop)

    # ------------------------------------------------------------
    # 13. Optionally write final MatrixTable
    # ------------------------------------------------------------
    if args.out_mt_path:
        mt_filtered.write(
            args.out_mt_path,
            overwrite=args.overwrite,
        )

    # ------------------------------------------------------------
    # 14. Export final VCF
    # ------------------------------------------------------------
    if args.metadata_vcf_or_header:
        metadata = hl.get_vcf_metadata(args.metadata_vcf_or_header)
    else:
        metadata = build_default_metadata()

    hl.export_vcf(
        mt_filtered,
        args.out_vcf,
        metadata=metadata
    )


if __name__ == "__main__":
    main()