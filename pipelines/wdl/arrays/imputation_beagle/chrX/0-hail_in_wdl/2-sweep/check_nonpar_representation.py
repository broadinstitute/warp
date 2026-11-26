#!/usr/bin/env python3
"""Show how a VDS represents genotypes in chrX non-PAR, by sex.

Reads one short window. Three parts:
  1. schema    entry / row / global fields of both halves of the VDS
  2. sparse    what is stored: ploidy and genotype class of variant-data calls, the
               call field of the reference data (if it has one), depth and GQ
  3. dense     what Hail produces when the window is densified (what Steps 1, 2
               and 4 of the pipeline see), on a subset of samples

Run on the Hail cluster once per cohort, a few minutes each:

    python3 check_nonpar_representation.py \\
        --vds gs://prod-drc-broad/v8/wgs/vds/aou_srwgs_short_variants_v8r1.vds/merged \\
        --sex-tsv gs://<bucket>/aou.foxtrot.sample_sex.tsv \\
        --temp-dir gs://<bucket>/tmp/repr_check

    python3 check_nonpar_representation.py \\
        --vds gs://fc-secure-759ddb91-68ac-4234-8fab-9c9325c3484f/CCDG_consent.vds \\
        --sex-tsv gs://<bucket>/anvil.sample_sex.tsv \\
        --temp-dir gs://<bucket>/tmp/repr_check

--sex-tsv is a two-column (sample, sex) table; without it everything is reported
under sex "NA". --dense-samples 0 skips part 3.
"""
import argparse
import re


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vds", required=True)
    ap.add_argument("--sex-tsv")
    ap.add_argument("--interval", default="chrX:3200000-3220000", help="a short non-PAR window")
    ap.add_argument("--dense-samples", type=int, default=2000, help="samples used for the densified part (0 = skip)")
    ap.add_argument("--ref-pad-bp", type=int, default=1000000,
                    help="densified part only: how far before the window to look for reference blocks "
                         "when the VDS has no stored maximum reference-block length")
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
    print("Hail version: %s" % hl.version())

    text = args.interval.replace(",", "")
    contig, coords = text.split(":")
    start, end = (int(x) for x in coords.split("-"))
    interval = hl.parse_locus_interval(text, reference_genome="GRCh38")

    vds = hl.vds.read_vds(args.vds)
    fd = hl.vds.VariantDataset.ref_block_max_length_field
    has_max_len = fd in vds.reference_data.globals

    # ------------------------------------------------------------------ 1. schema
    print("\n=== 1. SCHEMA ===")
    for name, mt in (("variant_data", vds.variant_data), ("reference_data", vds.reference_data)):
        print("%s entry fields : %s" % (name, ", ".join("%s:%s" % (f, mt[f].dtype) for f in mt.entry)))
        print("%s row fields   : %s" % (name, ", ".join(mt.row)))
        print("%s global fields: %s" % (name, ", ".join(mt.globals) or "(none)"))
    print("stored maximum reference-block length: %s"
          % (hl.eval(vds.reference_data.index_globals()[fd]) if has_max_len else "none"))

    # ------------------------------------------------------------------ sex labels
    sex_ht = None
    if args.sex_tsv:
        rows = []
        with hl.hadoop_open(args.sex_tsv, "r") as fh:
            for line in fh:
                parts = re.split(r"[\t, ]+", line.strip())
                if len(parts) >= 2:
                    rows.append({"s": parts[0], "sex": parts[1].upper()})
        sex_ht = hl.Table.parallelize(rows, schema=hl.tstruct(s=hl.tstr, sex=hl.tstr), key="s").distinct()

    def with_sex(mt):
        if sex_ht is None:
            return mt.annotate_cols(sex="NA")
        return mt.annotate_cols(sex=hl.coalesce(sex_ht[mt.s].sex, "NA"))

    def show(title, counts):
        print("\n-- " + title)
        total = sum(counts.values()) or 1
        for key, n in sorted(counts.items(), key=lambda kv: tuple(str(v) for v in kv[0].values())):
            print("   " + "  ".join("%s=%-8s" % (k, v) for k, v in key.items()) + " %14d  %6.2f%%" % (n, 100.0 * n / total))

    def show_stats(title, by_sex):
        print("\n-- " + title)
        for sex in sorted(by_sex, key=str):
            st = by_sex[sex]
            print("   sex=%-4s n=%-12d mean=%s" % (sex, st.n, "NA" if st.n == 0 else "%.2f" % st.mean))

    def call_class(c):
        """ref / het / hom_var for diploid calls; ref / alt for haploid ones."""
        return (hl.case()
                .when(hl.is_missing(c), "missing")
                .when(c.ploidy == 1, hl.if_else(c.is_non_ref(), "alt", "ref"))
                .when(c.is_hom_ref(), "ref")
                .when(c.is_het(), "het")
                .default("hom_var"))

    # ------------------------------------------------------------------ 2. sparse
    print("\n=== 2. SPARSE (as stored), %s, all samples ===" % text)
    vd = with_sex(hl.filter_intervals(vds.variant_data, [interval]))
    print("samples by sex: %s" % dict(vd.aggregate_cols(hl.agg.counter(vd.sex))))
    print("variant rows in window: %d" % vd.count_rows())

    gt = hl.vds.lgt_to_gt(vd.LGT, vd.LA) if "LGT" in vd.entry else vd.GT
    var_aggs = dict(
        calls=hl.agg.counter(hl.struct(sex=vd.sex, ploidy=gt.ploidy, cls=call_class(gt), phased=gt.phased)),
        gq=hl.agg.group_by(vd.sex, hl.agg.stats(vd.GQ)))
    if "LAD" in vd.entry:
        var_aggs["depth"] = hl.agg.group_by(vd.sex, hl.agg.stats(hl.sum(vd.LAD)))
    res = vd.aggregate_entries(hl.struct(**var_aggs))
    show("variant data: calls by sex, ploidy and class", res.calls)
    show_stats("variant data: GQ", res.gq)
    if "depth" in res:
        show_stats("variant data: depth at variant calls (sum of LAD)", res.depth)

    # Reference rows that START in the window are enough to see the representation
    rd = with_sex(hl.filter_intervals(vds.reference_data, [interval]))
    ref_call = "GT" if "GT" in rd.entry else ("LGT" if "LGT" in rd.entry else None)
    print("\nreference data call field: %s" % (ref_call or "none (densify fills diploid 0/0)"))
    block_len = rd.LEN if "LEN" in rd.entry else rd.END - rd.locus.position + 1
    ref_aggs = dict(length=hl.agg.group_by(rd.sex, hl.agg.stats(block_len)))
    if "GQ" in rd.entry:
        ref_aggs["gq"] = hl.agg.group_by(rd.sex, hl.agg.stats(rd.GQ))
    if "DP" in rd.entry:
        ref_aggs["dp"] = hl.agg.group_by(rd.sex, hl.agg.stats(rd.DP))
    if ref_call:
        ref_aggs["calls"] = hl.agg.counter(hl.struct(sex=rd.sex, ploidy=rd[ref_call].ploidy, cls=call_class(rd[ref_call])))
    res = rd.aggregate_entries(hl.struct(**ref_aggs))
    if "calls" in res:
        show("reference blocks: stored call by sex and ploidy", res.calls)
    show_stats("reference blocks: length (bp)", res.length)
    if "gq" in res:
        show_stats("reference blocks: GQ", res.gq)
    if "dp" in res:
        show_stats("reference blocks: DP", res.dp)

    # ------------------------------------------------------------------- 3. dense
    if args.dense_samples > 0:
        print("\n=== 3. DENSIFIED (what the pipeline's densify steps see), first %d samples ===" % args.dense_samples)
        pad = (int(hl.eval(vds.reference_data.index_globals()[fd])) - 1) if has_max_len else args.ref_pad_bp
        padded = hl.parse_locus_interval("%s:%d-%d" % (contig, max(1, start - pad), end), reference_genome="GRCh38")
        keep = vds.variant_data.cols().select().head(args.dense_samples)
        sub = hl.vds.filter_samples(
            hl.vds.VariantDataset(hl.filter_intervals(vds.reference_data, [padded]),
                                  hl.filter_intervals(vds.variant_data, [interval])),
            keep, keep=True)

        def dense_counts(ref_mt, label):
            mt = with_sex(hl.vds.to_dense_mt(hl.vds.VariantDataset(ref_mt, sub.variant_data)))
            c = hl.vds.lgt_to_gt(mt.LGT, mt.LA) if "LGT" in mt.entry else mt.GT
            show(label, mt.aggregate_entries(hl.agg.counter(hl.struct(sex=mt.sex, ploidy=c.ploidy, cls=call_class(c)))))

        dense_counts(sub.reference_data, "dense, VDS as stored (the original Step 1 / Step 2 scripts)")
        if ref_call:
            dense_counts(sub.reference_data.drop(ref_call),
                         "dense, reference call field dropped (the optimized Steps 1, 2 and 4)")


if __name__ == "__main__":
    main()
