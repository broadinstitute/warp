#!/usr/bin/env python3
"""Tabulate genotype ploidy right after densifying, the way Step 1 / Step 2 see it.

Answers: in non-PAR, are male reference genotypes diploid 0/0, haploid 0, or
missing once the VDS is densified? That decides whether Hail's excess-het test
in Step 1 is effectively females-only (haploid / missing) or diluted by male
homozygotes (diploid).

Run on the Hail cluster (small: one short window, a few thousand samples):

    python3 check_step1_ploidy.py \\
        --vds gs://prod-drc-broad/v8/wgs/vds/aou_srwgs_short_variants_v8r1.vds/merged \\
        --sex-tsv aou.foxtrot.sample_sex.tsv \\
        --interval chrX:3200000-3220000 \\
        --n-samples 2000 \\
        --temp-dir gs://<bucket>/tmp/ploidy_check

--sex-tsv (two columns: sample, sex) is optional; without it counts are not
split by sex. Run again with the AnVIL VDS to see its encoding.
"""
import argparse


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vds", required=True)
    ap.add_argument("--sex-tsv")
    ap.add_argument("--interval", default="chrX:3200000-3220000")
    ap.add_argument("--n-samples", type=int, default=2000)
    ap.add_argument("--ref-pad-bp", type=int, default=1000000)
    ap.add_argument("--temp-dir", required=True)
    ap.add_argument("--driver-memory", default="auto")
    args = ap.parse_args()

    import hail as hl
    driver_memory = args.driver_memory
    if driver_memory == "auto":
        with open("/proc/meminfo") as fh:
            driver_memory = "%dg" % max(1, int(int(fh.readline().split()[1]) / 1024 / 1024 * 0.7))
    hl.init(idempotent=True, tmp_dir=args.temp_dir.rstrip("/"), spark_conf={"spark.driver.memory": driver_memory})
    hl.default_reference("GRCh38")

    text = args.interval.replace(",", "")
    contig, coords = text.split(":")
    start, end = (int(x) for x in coords.split("-"))
    interval = hl.parse_locus_interval(text, reference_genome="GRCh38")
    padded = hl.parse_locus_interval("%s:%d-%d" % (contig, max(1, start - args.ref_pad_bp), end),
                                     reference_genome="GRCh38")

    vds = hl.vds.read_vds(args.vds)
    # Pre-filter both halves so a VDS without a stored maximum reference-block
    # length does not make Hail scan every reference partition.
    vds = hl.vds.VariantDataset(hl.filter_intervals(vds.reference_data, [padded]),
                                hl.filter_intervals(vds.variant_data, [interval]))
    keep = vds.variant_data.cols().select().head(args.n_samples)
    sample_ids = keep.s.collect()
    vds = hl.vds.filter_samples(vds, keep, keep=True)
    vds = hl.vds.filter_intervals(vds, [interval], split_reference_blocks=True)

    sex = {}
    if args.sex_tsv:
        wanted = set(sample_ids)
        opener = hl.hadoop_open if args.sex_tsv.startswith("gs://") else open
        with opener(args.sex_tsv, "r") as fh:
            for line in fh:
                f = line.split()
                if len(f) >= 2 and f[0] in wanted:
                    sex[f[0]] = f[1]
    print("samples: %d, with sex: %d" % (len(sample_ids), len(sex)))

    mt = hl.vds.to_dense_mt(vds)                      # same call as Step 1 / Step 2
    mt = mt.annotate_entries(GT=hl.vds.lgt_to_gt(mt.LGT, mt.LA))
    sex_lit = hl.literal(sex, hl.tdict(hl.tstr, hl.tstr))
    mt = mt.annotate_cols(sex=sex_lit.get(mt.s, "NA"))
    key = hl.struct(sex=mt.sex,
                    gt_defined=hl.is_defined(mt.GT),
                    ploidy=mt.GT.ploidy,
                    non_ref=mt.GT.is_non_ref())
    counts = mt.aggregate_entries(hl.agg.counter(key))
    n_rows = mt.count_rows()

    print("\nsites in window: %d" % n_rows)
    print("%-6s %-10s %-7s %-8s %12s" % ("sex", "defined", "ploidy", "non_ref", "entries"))
    for k, n in sorted(counts.items(), key=lambda kv: (str(kv[0].sex), str(kv[0].gt_defined), str(kv[0].ploidy), str(kv[0].non_ref))):
        print("%-6s %-10s %-7s %-8s %12d" % (k.sex, k.gt_defined, k.ploidy, k.non_ref, n))
    print("\nRead it as: XY rows with ploidy 2 and non_ref False are male reference genotypes "
          "that Hail's HWE test counts as homozygotes.")


if __name__ == "__main__":
    main()
