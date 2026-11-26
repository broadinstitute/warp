#!/usr/bin/env python3
"""Export the two run-specific site lists that vcfdist_report_sex.py needs:

    step1_pass.sites.tsv.bgz       sites in a Hail WDL run's Step 1 output
    anvil_exclusion.sites.tsv.bgz  sites on that run's Step 2 exclusion list

Nothing is densified and the VDSs are not read. Use the Step 1 and Step 2
outputs of the same Hail WDL run that produced the panel being evaluated, and
write into the folder that already holds aou.sites.tsv.bgz and
anvil.sites.tsv.bgz (those two do not depend on the run and are not touched).

Run on the Hail cluster:

    # outputs of the original pipeline (dense VCF and rows-only MT)
    python3 export_run_site_lists.py \\
        --step1-vcf gs://<run-output>/chrX.step1.aou.vcf.bgz \\
        --anvil-exclusion-mt gs://<run-output>/chrX.step2.anvil.mt \\
        ...
    # outputs of the optimized pipeline (Hail Tables)
    python3 export_run_site_lists.py \\
        --step1-ht gs://<run-output>/chrX.step1.aou.pass.ht \\
        --anvil-exclusion-ht gs://<run-output>/chrX.step2.anvil.exclusion.ht \\
        --interval chrX:1-5000000 \\
        --temp-dir gs://<bucket>/tmp/run_site_lists \\
        --out-dir  gs://<bucket>/site_tables_run2
"""
import argparse


def row_sites(hl, ht, interval):
    """contig / pos / ref / alts for a table keyed by locus and alleles."""
    ht = hl.filter_intervals(ht, [interval])
    ht = ht.select(contig=ht.locus.contig, pos=ht.locus.position,
                   ref=ht.alleles[0], alts=hl.delimit(ht.alleles[1:], ","))
    return ht.key_by().drop("locus", "alleles")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--step1-vcf", help="Step 1 output VCF (<output>/<contig>.step1.aou.vcf.bgz)")
    ap.add_argument("--anvil-exclusion-mt", help="Step 2 output MT (<output>/<contig>.step2.anvil.mt)")
    ap.add_argument("--step1-ht", help="optimized pipeline: Step 1 output table (<output>/<contig>.step1.aou.pass.ht)")
    ap.add_argument("--anvil-exclusion-ht", help="optimized pipeline: Step 2 output table (<output>/<contig>.step2.anvil.exclusion.ht)")
    ap.add_argument("--interval", default="chrX:1-5000000")
    ap.add_argument("--temp-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--driver-memory", default="auto", help="e.g. 60g; auto = 70%% of this machine's RAM")
    args = ap.parse_args()
    if not (args.step1_vcf or args.anvil_exclusion_mt or args.step1_ht or args.anvil_exclusion_ht):
        ap.error("give at least one of --step1-vcf, --step1-ht, --anvil-exclusion-mt, --anvil-exclusion-ht")
    out = args.out_dir.rstrip("/")

    import hail as hl
    driver_memory = args.driver_memory
    if driver_memory == "auto":
        with open("/proc/meminfo") as fh:
            driver_memory = "%dg" % max(1, int(int(fh.readline().split()[1]) / 1024 / 1024 * 0.7))
    hl.init(idempotent=True, tmp_dir=args.temp_dir.rstrip("/"), spark_conf={"spark.driver.memory": driver_memory})
    hl.default_reference("GRCh38")
    interval = hl.parse_locus_interval(args.interval.replace(",", ""), reference_genome="GRCh38")

    if args.step1_vcf:
        # Same import call as Step 3; only the row keys are used
        mt = hl.import_vcf(args.step1_vcf, reference_genome="GRCh38", force_bgz=True)
        row_sites(hl, mt.rows().select(), interval).export(out + "/step1_pass.sites.tsv.bgz")
        print("wrote " + out + "/step1_pass.sites.tsv.bgz")
    if args.anvil_exclusion_mt:
        mt = hl.read_matrix_table(args.anvil_exclusion_mt)
        row_sites(hl, mt.rows().select(), interval).export(out + "/anvil_exclusion.sites.tsv.bgz")
        print("wrote " + out + "/anvil_exclusion.sites.tsv.bgz")
    if args.step1_ht:
        # keyed by locus and alleles
        row_sites(hl, hl.read_table(args.step1_ht).select(), interval).export(out + "/step1_pass.sites.tsv.bgz")
        print("wrote " + out + "/step1_pass.sites.tsv.bgz")
    if args.anvil_exclusion_ht:
        # keyed by locus only; the report matches the exclusion list by position
        ht = hl.filter_intervals(hl.read_table(args.anvil_exclusion_ht), [interval])
        ht = ht.select(contig=ht.locus.contig, pos=ht.locus.position, ref="", alts="").key_by().drop("locus")
        ht.export(out + "/anvil_exclusion.sites.tsv.bgz")
        print("wrote " + out + "/anvil_exclusion.sites.tsv.bgz")


if __name__ == "__main__":
    main()
