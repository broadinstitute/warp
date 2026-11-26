"""Steps 1 and 2: per-site QC statistics for one cohort, for every site in the region.

Run once with --cohort aou (Step 1) and once with --cohort anvil (Step 2). Nothing
is filtered here: the thresholds are applied in Step 3, so trying a different
filter set does not need this step again, and filter_sweep.py can score filter
sets offline from the exported table.

Output: a Hail Table keyed by locus and alleles, and the same rows as a TSV.
    region                PAR1 | nonPAR | PAR2 | autosome
    n_alleles, max_ac     alleles in the row; largest alt allele count
    n_ad, mean_sum_ad     calls with AD defined; their mean total depth ("AVSAD")
    call_rate, mean_gq    over all samples, reference blocks included
    het_test_set          samples in the excess-het test: XX in chrX non-PAR, all elsewhere
    n_test_diploid        test-set samples with a diploid call at the site (densified
                          reference blocks included; missing and haploid calls are not)
    eh_n_het, eh_n_hom_var, eh_p, eh_f
                          one value per alt allele, in ALT order. Each alt allele is
                          tested against "every other allele" (what Hail's test gives
                          after splitting multiallelics): test-set diploid calls with
                          one copy / two copies of the allele, Hail's one-sided exact
                          excess-het p-value, and the inbreeding coefficient
                          F = 1 - observed het / expected het (negative = het excess)
    p_excess_het, excess_het_f
                          site summary: the smallest eh_p, and eh_f of that allele
    n_xy_called, n_xy_nonref, n_xy_het, n_xy_haploid_nonref, xy_het_frac
                          chrX only: XY samples with a call / a non-reference call /
                          a diploid het call / a haploid non-reference call, and
                          xy_het_frac = n_xy_het / n_xy_nonref. In non-PAR this is the
                          male het fraction: 0 if male calls are haploid, and otherwise
                          the share of male carriers called heterozygous.
    aou only:   f_LowQual, f_NO_HQ_GENOTYPES, f_ExcessHet   (precomputed site flags)
    anvil only: n_hq      samples passing the original NO_HQ_GENOTYPES condition

How the work is split:
  * a sparse pass over the variant data writes everything that depends only on
    variant calls (depth, allele counts, XY het counts, flags);
  * a dense pass adds what needs every sample (call rate, mean GQ, hom-ref counts,
    n_hq). Sites whose largest alt allele count is below --dense_min_ac are left
    out of the dense pass (their dense columns are missing): with the default of 2
    for AoU that is every site the allele-count filter removes anyway;
  * interval filtering is pushed down to both halves of the VDS by hand, so a VDS
    without a stored maximum reference-block length is not scanned in full.
"""
import argparse
import re
import subprocess
import sys
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument('--cohort',                     required=True, choices=['aou', 'anvil'])
parser.add_argument('--contig',                     required=True, help='contig to run, in hg38 format (e.g. chr20)')
parser.add_argument('--temp_bucket',                required=True, help='GCS path used for Hail temp files, likely the cluster temp bucket')
parser.add_argument('--bed_url',                    required=True, help='GCS url to BED file for interval filtering')
parser.add_argument('--sex_tsv_url',                default='',    help='GCS url to a two-column (sample, sex) table for this cohort, with sample IDs as they are in the VDS. Required for chrX')
parser.add_argument('--ref_pad_bp',                 type=int, default=1000000, help='how far before each interval to look for reference blocks when the VDS has no stored maximum reference-block length')
parser.add_argument('--dense_min_ac',               type=int, default=-1, help='sites with max alt AC below this skip the dense pass; default 2 for aou, 0 (every site) for anvil')
parser.add_argument('--executor_memory',            default='22g')
parser.add_argument('--executor_cores',             default='4')
parser.add_argument('--rig_cost',                   type=float, default=0.0, help='cluster cost per hour, for the report only')
# INPUTS
parser.add_argument('--input_vds_url',              required=True, help='GCS url to the cohort VDS')
# OUTPUTS
parser.add_argument('--output_stats_ht_url',        required=True, help='GCS url where the site statistics table is written')
parser.add_argument('--output_stats_tsv_url',       required=True, help='GCS url where the same table is exported as TSV (.tsv.bgz)')
parser.add_argument('--output_report_url',          required=True, help='GCS url where the report output will be copied')
args = parser.parse_args()

IS_AOU = args.cohort == 'aou'
DENSE_MIN_AC = args.dense_min_ac if args.dense_min_ac >= 0 else (2 if IS_AOU else 0)
ALT_MAX = 31
NONPAR_START, NONPAR_END = 2781480, 155701382      # GRCh38 chrX
SEX_LABELS = {'XX': 'XX', 'F': 'XX', 'FEMALE': 'XX', 'XY': 'XY', 'M': 'XY', 'MALE': 'XY'}


def run_subprocess(cmd, errorMessage):
    if isinstance(cmd, list):
        cmd = ' '.join(cmd)
    try:
        return subprocess.check_output(cmd, shell=True, universal_newlines=True)
    except subprocess.CalledProcessError as e:
        print(errorMessage)
        print("Exited with " + str(e.returncode) + "-" + e.output)
        exit(1)


start = datetime.now()

import hail as hl  # noqa: E402

step = 'Step_1' if IS_AOU else 'Step_2'
temp_dir = f'{args.temp_bucket}/{step}/temp'
hl.init(idempotent=True,
        tmp_dir=temp_dir,
        spark_conf={"spark.executor.memory": args.executor_memory,
                    "spark.executor.cores": args.executor_cores,
                    "spark.driver.memory": "60g",
                    "spark.driver.cores": "32",
                    "spark.local.dir": temp_dir})
hl.default_reference('GRCh38')

CHR_VAR = args.contig
IS_CHRX = CHR_VAR == 'chrX'
LABEL = 'AoU' if IS_AOU else 'AnVIL'


# ---------------------------------------------------------------- shared helpers
def load_intervals(bed_url, contig):
    """BED intervals on one contig, as Python hl.Interval objects sorted by start."""
    bed = hl.import_bed(bed_url, reference_genome='GRCh38')
    ivs = bed.filter(bed.interval.start.contig == contig).interval.collect()
    if not ivs:
        sys.exit(f"No intervals on {contig} in {bed_url}")
    return sorted(ivs, key=lambda iv: iv.start.position)


def read_region(vds_path, contig, intervals, ref_pad_bp, label):
    """Read a VDS restricted to the intervals, pruning partitions on both halves.

    Variant data: rows inside the intervals. Reference data: rows that start inside
    the intervals or up to `pad` bp before them, so that blocks which begin earlier
    and reach into an interval are still there for densifying. `pad` is the stored
    maximum reference-block length when the VDS has one, otherwise ref_pad_bp.
    hl.vds.filter_intervals is deliberately not used: without a stored maximum it
    leaves the whole contig's reference data in the plan (a full pass).
    """
    vds = hl.vds.read_vds(vds_path)
    rd, vd = vds.reference_data, vds.variant_data
    fd = hl.vds.VariantDataset.ref_block_max_length_field
    if fd in rd.globals:
        pad = int(hl.eval(rd.index_globals()[fd])) - 1
        print(f"{label}: stored maximum reference-block length -> padding {pad} bp")
    else:
        pad = ref_pad_bp
        print(f"{label}: no stored maximum reference-block length -> padding {pad} bp (--ref_pad_bp)")
    padded = [hl.Interval(hl.Locus(contig, max(1, iv.start.position - pad), reference_genome='GRCh38'),
                          iv.end, includes_start=True, includes_end=iv.includes_end)
              for iv in intervals]
    rd = hl.filter_intervals(rd, padded)
    vd = hl.filter_intervals(vd, intervals)
    try:
        print(f"{label}: {vd.n_partitions()} variant partitions, {rd.n_partitions()} reference partitions after interval filter")
    except Exception as exc:  # informational only
        print(f"{label}: could not count partitions ({exc})")
    return rd, vd


def ref_entry_fields(rd, wanted):
    """Reference entry fields to keep: block end/length (whichever this Hail version stores) plus `wanted`."""
    return [f for f in ('END', 'LEN') if f in rd.entry] + [f for f in wanted if f in rd.entry]


def load_sex_table(sex_tsv_url):
    """Keyed table sample -> 'XX' | 'XY' from a two-column (sample, sex) file; other labels are skipped."""
    rows, n_rows = {}, 0
    with hl.hadoop_open(sex_tsv_url, 'r') as fh:
        for line in fh:
            parts = re.split(r'[\t, ]+', line.strip())
            if len(parts) < 2:
                continue
            n_rows += 1
            sex = SEX_LABELS.get(parts[1].upper())
            if sex:
                rows[parts[0]] = sex
    print(f"{sex_tsv_url}: {n_rows} rows, {sum(v == 'XX' for v in rows.values())} XX, {sum(v == 'XY' for v in rows.values())} XY")
    return hl.Table.parallelize([{'s': s, 'sex': x} for s, x in sorted(rows.items())],
                                schema=hl.tstruct(s=hl.tstr, sex=hl.tstr), key='s')


# NOTE: a Hail expression belongs to the table object it was built from, so the
# helpers below take the current table and are called again after every step that
# returns a new one.
def region_label(locus):
    if not IS_CHRX:
        return hl.str('autosome')
    return (hl.case()
            .when(locus.position < NONPAR_START, 'PAR1')
            .when(locus.position > NONPAR_END, 'PAR2')
            .default('nonPAR'))


def in_het_test(t):
    """Samples counted in the excess-het test at this row: XX in chrX non-PAR, everyone elsewhere."""
    if not IS_CHRX:
        return hl.literal(True)
    return hl.if_else((t.locus.position >= NONPAR_START) & (t.locus.position <= NONPAR_END), t.is_xx, True)


def global_gt(t):
    return hl.vds.lgt_to_gt(t.LGT, t.LA)


def sum_ad(t):
    # Sum of the global AD array == sum of the local one (local_to_global only adds zeros)
    return hl.or_missing(hl.is_defined(t.LA), hl.sum(t.LAD))


def hq_variant_entry(t):
    """The original Step 2 NO_HQ_GENOTYPES_COND for an entry of the variant data (missing -> False)."""
    gt = global_gt(t)
    ad = hl.vds.local_to_global(t.LAD, t.LA, n_alleles=hl.len(t.alleles), fill_value=0, number='R')
    s = hl.sum(ad)
    return hl.coalesce(
        (t.GQ >= 20) & (t.DP >= 10) &
        ((gt.is_het_ref() & ((ad[1] / s) >= 0.2)) |
         (gt.is_het_non_ref() & ((ad[0] / s) >= 0.2) & ((ad[1] / s) >= 0.2)) |
         (~gt.is_het_ref() & ~gt.is_het_non_ref())),
        False)


# ------------------------------------------------------------------------- read
intervals = load_intervals(args.bed_url, CHR_VAR)
rd, vd = read_region(args.input_vds_url, CHR_VAR, intervals, args.ref_pad_bp, LABEL)
n_samples = vd.count_cols()

if IS_CHRX:
    if not args.sex_tsv_url:
        sys.exit("--sex_tsv_url is required for chrX")
    sex_ht = load_sex_table(args.sex_tsv_url)
    vd = vd.annotate_cols(_sex=sex_ht[vd.s].sex)
    vd = vd.annotate_cols(is_xx=hl.coalesce(vd._sex == 'XX', False), is_xy=hl.coalesce(vd._sex == 'XY', False))
    n_xx, n_xy = vd.aggregate_cols((hl.agg.count_where(vd.is_xx), hl.agg.count_where(vd.is_xy)))
    print(f"{LABEL}: {n_samples} samples, {n_xx} XX, {n_xy} XY, {n_samples - n_xx - n_xy} without a usable label")
    if n_xx == 0 or n_xy == 0:
        sys.exit("No VDS sample matched an XX (or an XY) label in --sex_tsv_url; check the sample ID column")

# Sites with more than ALT_MAX alt alleles are never used (same cut as the original Steps 1, 2 and 5b)
vd = vd.filter_rows(hl.len(vd.alleles) < (ALT_MAX + 2))

# ------------------------------------------------------------- 1. sparse pass
gt = global_gt(vd)
depth = sum_ad(vd)
cs_all = hl.agg.call_stats(gt, vd.alleles)
# Copies of each alt allele among the test set's diploid calls. Reference blocks carry
# no alt allele, so these counts are complete without densifying.
copies = gt.one_hot_alleles(hl.len(vd.alleles))
per_alt = hl.agg.filter(
    in_het_test(vd) & (gt.ploidy == 2),
    hl.agg.array_agg(
        lambda k: hl.struct(n1=hl.agg.count_where(copies[k] == 1), n2=hl.agg.count_where(copies[k] == 2)),
        hl.range(1, hl.len(vd.alleles))))
sparse_fields = dict(
    region=region_label(vd.locus),
    n_alleles=hl.len(vd.alleles),
    max_ac=hl.max(cs_all.AC[1:]),
    n_ad=hl.agg.count_where(depth >= 0),
    mean_sum_ad=hl.agg.sum(depth) / hl.agg.count_where(depth >= 0),
    _per_alt=per_alt)
if IS_CHRX:
    sparse_fields.update(
        n_xy_nonref=hl.agg.count_where(vd.is_xy & gt.is_non_ref()),
        n_xy_het=hl.agg.count_where(vd.is_xy & (gt.ploidy == 2) & gt.is_het()),
        n_xy_haploid_nonref=hl.agg.count_where(vd.is_xy & (gt.ploidy == 1) & gt.is_non_ref()))
if IS_AOU:
    sparse_fields.update(
        f_LowQual=vd.filters.contains('LowQual'),
        f_NO_HQ_GENOTYPES=vd.filters.contains('NO_HQ_GENOTYPES'),
        f_ExcessHet=vd.filters.contains('ExcessHet'))
sparse_path = f'{temp_dir}/site_stats_sparse.ht'
vd.select_rows(**sparse_fields).rows().write(sparse_path, overwrite=True)
sp = hl.read_table(sparse_path)
print(f"{LABEL}: {sp.count()} sites; sparse statistics written")

# -------------------------------------------------------------- 2. dense pass
if DENSE_MIN_AC > 0:
    vd = vd.filter_rows(hl.coalesce(sp[vd.locus, vd.alleles].max_ac >= DENSE_MIN_AC, False))
# LGT is kept only where the global GT is defined, so "called" and "diploid hom-ref"
# mean the same as they did on GT in the original scripts.
dense_entries = dict(LGT=hl.or_missing(hl.is_defined(global_gt(vd)), vd.LGT), GQ=vd.GQ)
if not IS_AOU:
    dense_entries.update(DP=vd.DP, _hq=hq_variant_entry(vd))
vd = vd.select_entries(**dense_entries)
# No call field is kept in the reference data, so densify fills reference blocks with diploid 0/0
rd = rd.select_entries(*ref_entry_fields(rd, ['GQ'] if IS_AOU else ['GQ', 'DP']))
mt = hl.vds.to_dense_mt(hl.vds.VariantDataset(rd, vd))

dense_fields = dict(
    _n_called=hl.agg.count_where(hl.is_defined(mt.LGT)),
    mean_gq=hl.agg.stats(mt.GQ).mean,
    n_test_diploid=hl.agg.count_where(in_het_test(mt) & (mt.LGT.ploidy == 2)))
if IS_CHRX:
    dense_fields.update(n_xy_called=hl.agg.count_where(mt.is_xy & hl.is_defined(mt.LGT)))
if not IS_AOU:
    # For a densified reference entry (0/0, no AD) the original condition reduces to GQ >= 20 & DP >= 10
    dense_fields.update(n_hq=hl.agg.count_where(hl.coalesce(mt._hq, (mt.GQ >= 20) & (mt.DP >= 10))))
dense_path = f'{temp_dir}/site_stats_dense.ht'
mt.select_rows(**dense_fields).rows().write(dense_path, overwrite=True)
de = hl.read_table(dense_path)
print(f"{LABEL}: dense statistics written for {de.count()} sites")

# ------------------------------------------------------ 3. join and finish
ht = sp.annotate(**de[sp.key])

# Per alt allele: hom-other / het / hom counts among the test set's diploid calls, then
# Hail's one-sided exact test. With no test-set carrier the aggregation above is missing,
# which means zero carriers of every alt allele.
n_alt = ht.n_alleles - 1
counts = hl.coalesce(ht._per_alt, hl.range(n_alt).map(lambda _: hl.struct(n1=hl.int64(0), n2=hl.int64(0))))


def allele_test(c):
    n0 = hl.max(ht.n_test_diploid - c.n1 - c.n2, 0)
    hwe = hl.hardy_weinberg_test(hl.int32(n0), hl.int32(c.n1), hl.int32(c.n2), one_sided=True)
    informative = (c.n1 + c.n2 > 0) & (n0 + c.n1 > 0)       # the allele is seen, and so is another one
    het_obs = c.n1 / ht.n_test_diploid
    return hl.struct(
        n_het=c.n1,
        n_hom_var=c.n2,
        p=hl.or_missing(informative, hwe.p_value),
        f=hl.or_missing(informative & (hwe.het_freq_hwe > 0), 1 - het_obs / hwe.het_freq_hwe))


tests = hl.or_missing(hl.is_defined(ht.n_test_diploid), counts.map(allele_test))
ht = ht.annotate(_tests=tests)
ranked = hl.sorted(ht._tests.filter(lambda t: hl.is_defined(t.p)), key=lambda t: t.p)
ht = ht.annotate(
    call_rate=hl.float(ht._n_called) / n_samples,
    het_test_set=hl.if_else(ht.region == 'nonPAR', 'XX', 'all'),
    eh_n_het=ht._tests.map(lambda t: t.n_het),
    eh_n_hom_var=ht._tests.map(lambda t: t.n_hom_var),
    eh_p=ht._tests.map(lambda t: t.p),
    eh_f=ht._tests.map(lambda t: t.f),
    p_excess_het=hl.or_missing(hl.len(ranked) > 0, ranked[0].p),
    excess_het_f=hl.or_missing(hl.len(ranked) > 0, ranked[0].f))
if IS_CHRX:
    ht = ht.annotate(xy_het_frac=hl.or_missing(ht.n_xy_nonref > 0, ht.n_xy_het / ht.n_xy_nonref))
ht = ht.drop('_n_called', '_per_alt', '_tests')
ht.write(args.output_stats_ht_url, overwrite=True)

ht = hl.read_table(args.output_stats_ht_url)
flat = ht.key_by()
flat = flat.select(contig=flat.locus.contig, pos=flat.locus.position, ref=flat.alleles[0],
                   alts=hl.delimit(flat.alleles[1:], ','),
                   **{f: (hl.delimit(flat[f].map(lambda x: hl.coalesce(hl.str(x), 'NA')), ',')
                          if isinstance(flat[f].dtype, hl.tarray) else flat[f])
                      for f in flat.row if f not in ('locus', 'alleles')})
flat.export(args.output_stats_tsv_url)

# ----------------------------------------------------------------------- report
end = datetime.now()
time_est = (end - start).total_seconds() / 3600
export_report = (
    f"Chro:\t{CHR_VAR}\n"
    f"Path:\t{args.output_stats_ht_url}\n"
    f"Time:\t{time_est}\n"
    f"Rigs:\t{args.rig_cost}\n"
    f"Cost:\t{args.rig_cost * time_est}")
local_fname = "tmp_step_report.txt"
with open(local_fname, 'w') as fp:
    fp.write(export_report)
run_subprocess(f"gcloud storage cp {local_fname} {args.output_report_url}",
               f"Error copying report - {local_fname} to {args.output_report_url}")
