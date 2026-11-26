#!/usr/bin/env python3
"""Infer each sample's X-chromosome copy number class (XX or XY) from a VDS, and write
a two-column sex table that the pipeline accepts as aou_sex_tsv_url / anvil_sex_tsv_url.

Only sparse variant data is read (no reference blocks, no densifying), from many short
windows spread along chrX non-PAR plus a few on an autosome. Two independent signals
are computed per sample:

    depth ratio   mean depth at non-reference calls on chrX non-PAR divided by the same
                  on the autosome: about 1 with two X copies, about 0.5 with one
    het fraction  share of chrX non-PAR non-reference calls that are diploid het:
                  high with two X copies, near 0 with one (or all calls are haploid)

Windows are spread out on purpose. Within one contiguous stretch a woman's two X
haplotypes are often alike, so her het fraction there says little; over a hundred
separate loci it is stable.

    XY         depth ratio below --depth-split, and calls mostly haploid or het
               fraction at most --xy-max-het
    XX         depth ratio at or above --depth-split, and het fraction at least --xx-min-het
    ambiguous  the two signals disagree, or the het fraction is in between
    unknown    too few calls on chrX or on the autosome

"XX" means "behaves as two X copies in non-PAR", which is what the excess-het test
needs; it is not a karyotype. The pipeline uses XX and XY and ignores the rest.

    python3 infer_sex_from_vds.py \\
        --vds gs://fc-secure-759ddb91-68ac-4234-8fab-9c9325c3484f/CCDG_consent.vds \\
        --out gs://<bucket>/resources/anvil.inferred_sex.tsv \\
        --temp-dir gs://<bucket>/tmp/infer_sex

To check the method, run it on the AoU VDS with --compare-sex-tsv <the AoU sex table>:
it prints inferred class against the given label.
"""
import argparse
import collections
import re

NONPAR_START, NONPAR_END = 2781480, 155701382      # GRCh38 chrX
CONTIG_LENGTH = {"chr20": 64444167, "chr21": 46709983, "chr22": 50818468, "chr19": 58617616, "chr1": 248956422}


def spread(contig, first, last, n, width):
    """n windows of `width` bp, evenly spaced between first and last."""
    if n <= 0:
        return []
    step = (last - first - width) / max(1, n - 1)
    return ["%s:%d-%d" % (contig, int(first + i * step), int(first + i * step) + width - 1) for i in range(n)]


def histogram(values, title, lo=0.0, hi=1.0, bins=20):
    print("\n" + title)
    counts = collections.Counter()
    for v in values:
        counts[min(bins - 1, max(0, int((v - lo) / (hi - lo) * bins)))] += 1
    top = max(counts.values()) if counts else 1
    for b in range(bins):
        a = lo + b * (hi - lo) / bins
        print("  %.2f-%.2f %9d %s" % (a, a + (hi - lo) / bins, counts[b], "#" * int(50 * counts[b] / top)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vds", required=True)
    ap.add_argument("--out", required=True, help="gs:// path of the table to write (tab-separated, with a header)")
    ap.add_argument("--x-windows", type=int, default=100, help="windows spread along chrX non-PAR")
    ap.add_argument("--autosome-windows", type=int, default=20, help="windows spread along the autosome")
    ap.add_argument("--autosome", default="chr20", choices=sorted(CONTIG_LENGTH))
    ap.add_argument("--window-bp", type=int, default=10000)
    ap.add_argument("--depth-split", type=float, default=0.75, help="X / autosome depth ratio separating one X copy from two")
    ap.add_argument("--xx-min-het", type=float, default=0.35)
    ap.add_argument("--xy-max-het", type=float, default=0.15)
    ap.add_argument("--min-calls", type=int, default=50, help="minimum non-reference calls on chrX (and a fifth of that on the autosome)")
    ap.add_argument("--compare-sex-tsv", help="optional two-column (sample, sex) table to cross-tabulate against")
    ap.add_argument("--temp-dir", required=True)
    ap.add_argument("--driver-memory", default="auto")
    args = ap.parse_args()

    windows = (spread("chrX", NONPAR_START + 20000, NONPAR_END - 20000, args.x_windows, args.window_bp)
               + spread(args.autosome, 1000000, CONTIG_LENGTH[args.autosome] - 1000000, args.autosome_windows, args.window_bp))
    print("%d chrX non-PAR windows and %d %s windows of %d bp" % (args.x_windows, args.autosome_windows, args.autosome, args.window_bp))

    import hail as hl
    driver_memory = args.driver_memory
    if driver_memory == "auto":
        with open("/proc/meminfo") as fh:
            driver_memory = "%dg" % max(1, int(int(fh.readline().split()[1]) / 1024 / 1024 * 0.7))
    hl.init(idempotent=True, tmp_dir=args.temp_dir.rstrip("/"), spark_conf={"spark.driver.memory": driver_memory})
    hl.default_reference("GRCh38")

    intervals = [hl.parse_locus_interval(w, reference_genome="GRCh38") for w in windows]
    vd = hl.filter_intervals(hl.vds.read_vds(args.vds).variant_data, intervals)
    gt = hl.vds.lgt_to_gt(vd.LGT, vd.LA) if "LGT" in vd.entry else vd.GT
    if "LAD" in vd.entry:
        depth = hl.sum(vd.LAD)
    elif "AD" in vd.entry:
        depth = hl.sum(vd.AD)
    else:
        depth = vd.DP
    on_x = vd.locus.contig == "chrX"
    nonref = gt.is_non_ref()
    vd = vd.annotate_cols(
        x_nonref=hl.agg.count_where(on_x & nonref),
        x_het=hl.agg.count_where(on_x & (gt.ploidy == 2) & gt.is_het()),
        x_haploid=hl.agg.count_where(on_x & (gt.ploidy == 1) & nonref),
        x_depth=hl.agg.filter(on_x & nonref, hl.agg.mean(depth)),
        a_nonref=hl.agg.count_where(~on_x & nonref),
        a_depth=hl.agg.filter(~on_x & nonref, hl.agg.mean(depth)))
    ht = vd.cols().key_by()
    ht = ht.select(sample_id=ht.s, x_nonref=ht.x_nonref, x_het=ht.x_het, x_haploid=ht.x_haploid,
                   x_depth=ht.x_depth, a_nonref=ht.a_nonref, a_depth=ht.a_depth)
    ht = ht.checkpoint(args.temp_dir.rstrip("/") + "/infer_sex_counts.ht", overwrite=True)

    het_frac = hl.or_missing(ht.x_nonref > 0, ht.x_het / ht.x_nonref)
    ratio = hl.or_missing(ht.a_depth > 0, ht.x_depth / ht.a_depth)
    one_copy = ratio < args.depth_split
    sex = (hl.case()
           .when((ht.x_nonref < args.min_calls) | (ht.a_nonref < args.min_calls / 5) | hl.is_missing(ratio), "unknown")
           .when(one_copy & ((ht.x_haploid > 0.5 * ht.x_nonref) | (het_frac <= args.xy_max_het)), "XY")
           .when(~one_copy & (het_frac >= args.xx_min_het), "XX")
           .default("ambiguous"))
    ht = ht.select(sample_id=ht.sample_id, sex=sex, depth_ratio=ratio, het_frac=het_frac,
                   x_nonref=ht.x_nonref, x_haploid=ht.x_haploid, a_nonref=ht.a_nonref)
    ht.export(args.out)

    rows = [r for r in ht.collect()]
    known = [r for r in rows if r.sex != "unknown"]
    print("\n%d samples" % len(rows))
    print("inferred class: %s" % dict(collections.Counter(r.sex for r in rows)))
    histogram([r.depth_ratio for r in known], "chrX non-PAR / %s depth ratio (two X copies near 1.0, one near 0.5):" % args.autosome, 0.0, 2.0)
    histogram([r.het_frac for r in known if r.depth_ratio >= args.depth_split],
              "het fraction, samples with depth ratio >= %.2f (two copies; expect one high peak):" % args.depth_split)
    histogram([r.het_frac for r in known if r.depth_ratio < args.depth_split],
              "het fraction, samples with depth ratio < %.2f (one copy; expect a peak at 0):" % args.depth_split)

    if args.compare_sex_tsv:
        given = {}
        with hl.hadoop_open(args.compare_sex_tsv, "r") as fh:
            for line in fh:
                parts = re.split(r"[\t, ]+", line.strip())
                if len(parts) >= 2:
                    given[parts[0]] = parts[1].upper()
        cross = collections.Counter((given.get(r.sample_id, "not listed"), r.sex) for r in rows)
        print("\ngiven label x inferred class:")
        for (g, i), n in sorted(cross.items()):
            print("  %-12s %-10s %9d" % (g, i, n))

    print("\nWrote " + args.out)


if __name__ == "__main__":
    main()
