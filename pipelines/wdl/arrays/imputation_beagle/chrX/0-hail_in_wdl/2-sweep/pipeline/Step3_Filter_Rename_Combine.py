"""Step 3: apply the site filters, restrict both VDSs to the kept sites, rename AnVIL
samples, and combine.

  * The filters are applied here, to the site-statistics tables from Steps 1 and 2,
    so a different filter set only needs this step and Step 4 again. --filter_spec
    sets the thresholds (see the block below); with no spec they are the original
    ones: AVSAD >= 12, call rate >= 0.90, mean GQ >= 30, AoU max alt AC >= 2, AoU
    LowQual / NO_HQ_GENOTYPES flags, AnVIL "no high-quality genotype", and Hail's
    excess-het test at p < 3.4e-6 (per alt allele, on diploid calls; XX samples in chrX
    non-PAR, all samples elsewhere).
  * kept sites = AoU sites that pass, minus loci on the AnVIL exclusion list (matched
    by locus, as before). The pass list, the exclusion list and the kept list are
    written out as small tables.
  * Each VDS is cut down to the kept sites and to the shared entry schema, AnVIL
    samples are renamed, and the two are merged with the same hl.vds.new_combiner
    call as the original Step 4.
  * This replaces the original Steps 3, 3b, 4 and 5a. Step 5a is dropped because
    both inputs to the combiner already contain only kept sites.
  * Interval filtering is pushed down to both halves of each VDS by hand. Reference
    blocks are not split at interval boundaries (densifying only looks at variant
    sites, all inside the intervals), but reference rows before an interval are kept
    only where a block reaches into it.
"""
import argparse
import subprocess
import sys
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument('--contig',                         required=True, help='contig to run, in hg38 format (e.g. chr20)')
parser.add_argument('--temp_bucket',                    required=True, help='GCS path used for Hail temp files, likely the cluster temp bucket')
parser.add_argument('--bed_url',                        required=True, help='GCS url to BED file for interval filtering')
parser.add_argument('--filter_spec',                    default='none', help='filter thresholds, e.g. aou.par.min_call_rate=0.8,anvil.nonpar.min_avsad=6 ; "none" = the defaults')
parser.add_argument('--ref_pad_bp',                     type=int, default=1000000, help='how far before each interval to look for reference blocks when a VDS has no stored maximum reference-block length')
parser.add_argument('--executor_memory',                default='22g')
parser.add_argument('--executor_cores',                 default='4')
parser.add_argument('--rig_cost',                       type=float, default=0.0, help='cluster cost per hour, for the report only')
# INPUTS
parser.add_argument('--input_aou_stats_ht_url',         required=True, help='AoU site statistics, from Step 1')
parser.add_argument('--input_anvil_stats_ht_url',       required=True, help='AnVIL site statistics, from Step 2')
parser.add_argument('--input_aou_vds_url',              required=True, help='GCS url to the AoU VDS')
parser.add_argument('--input_anvil_vds_url',            required=True, help='GCS url to the AnVIL CCDG VDS')
# OUTPUTS
parser.add_argument('--output_aou_pass_ht_url',         required=True, help='AoU sites that pass the filters')
parser.add_argument('--output_anvil_exclusion_ht_url',  required=True, help='AnVIL loci that fail the filters')
parser.add_argument('--output_keep_ht_url',             required=True, help='kept loci (pass list minus exclusion list)')
parser.add_argument('--output_aou_filtered_vds_url',    required=True, help='GCS url where the filtered AoU VDS is written')
parser.add_argument('--output_anvil_filtered_vds_url',  required=True, help='GCS url where the filtered, renamed AnVIL VDS is written')
parser.add_argument('--output_combined_vds_url',        required=True, help='GCS url where the combined AoU-AnVIL VDS is written')
parser.add_argument('--output_report_url',              required=True, help='GCS url where the report output will be copied')
args = parser.parse_args()

# Fields dropped / kept so that both cohorts share one schema (unchanged from the original Step 3)
DROP_AOU_VAR = ['truth_sensitivity_snp_threshold', 'truth_sensitivity_indel_threshold', 'filters', 'as_vets', 'PS', 'FT']
DROP_AOU_REF = ['ref_block_max_length', 'GT']
DROP_ANVIL_VAR = ['LPGT', 'LPL', 'gvcf_info', 'DP', 'MIN_DP', 'PID', 'SB', 'rsid']
DROP_ANVIL_REF = ['ref_allele', 'DP']
SELECT_VAR = ['GQ', 'RGQ', 'LGT', 'LAD', 'LA']
SELECT_REF = ['GQ', 'END']


# ----------------------------------------------------------------------------
# Filter specification (plain Python; the same block is in Step3_Filter_Rename_Combine.py
# and filter_sweep.py, so a spec that wins the sweep can be pasted into the WDL).
#
# A spec is a comma-separated list of key=value with no spaces, e.g.
#     aou.par.min_call_rate=0.8,anvil.nonpar.min_avsad=6,aou.excess_het=flag
# Keys are <cohort>.<name> or <cohort>.<region>.<name>
#     cohort  aou | anvil | both
#     region  par | nonpar          (omit to set the value everywhere, autosomes included)
#     name    min_avsad | min_call_rate | min_gq | min_ac (aou only) | excess_het_p |
#             excess_het_f (off by default: if set, the excess-het test only removes a site
#             whose inbreeding coefficient is also below this, e.g. -0.3) |
#             max_xy_het_frac (chrX non-PAR only; off by default: remove a site where more
#             than this share of XY carriers are called het) |
#             xy_het_min_carriers (max_xy_het_frac only applies to sites with at least
#             this many XY carriers; default 20)
# plus  aou.excess_het=hail|flag   hail: Hail's excess-het test everywhere, per alt allele, on
#                                        diploid calls (XX samples in non-PAR, all elsewhere)
#                                  flag: the same test in non-PAR, AoU's precomputed
#                                        ExcessHet site flag elsewhere
# A value of "none" switches a threshold off. An empty spec, or "none", gives the defaults.
FILTER_DEFAULTS = {
    'aou':   {'min_avsad': 12.0, 'min_call_rate': 0.9, 'min_gq': 30.0, 'min_ac': 2.0,
              'excess_het_p': 3.4e-6, 'excess_het_f': None, 'max_xy_het_frac': None,
              'xy_het_min_carriers': 20.0},
    'anvil': {'min_avsad': 12.0, 'min_call_rate': 0.9, 'min_gq': 30.0,
              'excess_het_p': 3.4e-6, 'excess_het_f': None, 'max_xy_het_frac': None,
              'xy_het_min_carriers': 20.0},
}
FILTER_REGIONS = ('all', 'par', 'nonpar')


def parse_filter_spec(spec):
    """-> ({cohort: {name: {'all': v, 'par': v, 'nonpar': v}}}, aou_excess_het_mode)."""
    thresholds = {c: {n: {r: v for r in FILTER_REGIONS} for n, v in d.items()}
                  for c, d in FILTER_DEFAULTS.items()}
    mode = 'hail'
    spec = (spec or '').strip()
    if spec.lower() in ('', 'none', 'default'):
        return thresholds, mode
    items = []
    for item in spec.split(','):
        if '=' not in item:
            raise ValueError("filter spec item %r is not key=value" % item)
        key, value = item.split('=', 1)
        items.append((key.strip().lower().split('.'), value.strip().lower(), item))
    # general settings first, region-specific ones last, so the more specific key wins
    for parts, value, item in sorted(items, key=lambda x: len(x[0])):
        if parts == ['aou', 'excess_het']:
            if value not in ('hail', 'flag'):
                raise ValueError("aou.excess_het must be hail or flag, not %r" % value)
            mode = value
            continue
        if len(parts) == 2:
            cohort, regions, name = parts[0], FILTER_REGIONS, parts[1]
        elif len(parts) == 3 and parts[1] in ('par', 'nonpar'):
            cohort, regions, name = parts[0], (parts[1],), parts[2]
        else:
            raise ValueError("cannot read filter spec key in %r" % item)
        cohorts = ('aou', 'anvil') if cohort == 'both' else (cohort,)
        number = None if value == 'none' else float(value)
        for c in cohorts:
            if c not in thresholds or name not in thresholds[c]:
                if cohort == 'both' and c in thresholds:
                    continue          # e.g. both.min_ac: AnVIL has no such filter
                raise ValueError("unknown filter %r in %r" % ('.'.join(parts), item))
            for r in regions:
                thresholds[c][name][r] = number
    return thresholds, mode
# ----------------------------------------------------------------------------

THRESHOLDS, AOU_EXCESS_HET = parse_filter_spec(args.filter_spec)   # fail on a bad spec before starting Hail
print(f"filter spec: {args.filter_spec}")
for _cohort in ('aou', 'anvil'):
    for _name, _by_region in sorted(THRESHOLDS[_cohort].items()):
        print(f"  {_cohort}.{_name}: {_by_region}")
print(f"  aou.excess_het: {AOU_EXCESS_HET}")


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

temp_dir = f'{args.temp_bucket}/Step_3/temp'
hl.init(idempotent=True,
        tmp_dir=temp_dir,
        spark_conf={"spark.executor.memory": args.executor_memory,
                    "spark.executor.cores": args.executor_cores,
                    "spark.driver.memory": "60g",
                    "spark.driver.cores": "32",
                    "spark.local.dir": temp_dir})
hl.default_reference('GRCh38')

CHR_VAR = args.contig


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


def drop_present(mt, names):
    """mt.drop for the names that exist (row, entry, column or global fields)."""
    present = [n for n in names if n in mt.row or n in mt.entry or n in mt.col or n in mt.globals]
    return mt.drop(*present) if present else mt


def block_end(rd):
    """Last position covered by a reference block, whichever field this Hail version stores."""
    return rd.END if 'END' in rd.entry else rd.locus.position + rd.LEN - 1


# ---------------------------------------------------------------- site filters
def threshold(ht, cohort, name):
    """Region-aware threshold as a Hail expression (missing where the threshold is off)."""
    by_region = THRESHOLDS[cohort][name]
    value = {r: (hl.missing(hl.tfloat64) if v is None else hl.float64(v)) for r, v in by_region.items()}
    return (hl.case()
            .when(ht.region == 'nonPAR', value['nonpar'])
            .when(ht.region.startswith('PAR'), value['par'])
            .default(value['all']))


def below(ht, x, cohort, name):
    """x < threshold; False where the threshold is switched off."""
    t = threshold(ht, cohort, name)
    return hl.if_else(hl.is_defined(t), x < t, False)


def above(ht, x, cohort, name):
    """x > threshold; False where the threshold is switched off."""
    t = threshold(ht, cohort, name)
    return hl.if_else(hl.is_defined(t), x > t, False)


def common_clauses(ht, cohort):
    """Clauses shared by both cohorts. Each is True where the site fails; the low-depth,
    call-rate and GQ clauses are written as in the original scripts."""
    clauses = [
        (ht.n_ad >= 1) & below(ht, ht.mean_sum_ad, cohort, 'min_avsad'),
        below(ht, ht.call_rate, cohort, 'min_call_rate'),
        below(ht, ht.mean_gq, cohort, 'min_gq')]
    if 'xy_het_frac' in ht.row:
        # male het fraction, chrX non-PAR only; off unless max_xy_het_frac is set. Sites with
        # fewer than xy_het_min_carriers XY carriers are left alone: one stray het call in a
        # single male gives a fraction of 1 and says nothing about the site.
        min_carriers = threshold(ht, cohort, 'xy_het_min_carriers')
        enough = hl.if_else(hl.is_defined(min_carriers), hl.coalesce(ht.n_xy_nonref >= min_carriers, False), True)
        clauses.append((ht.region == 'nonPAR') & enough &
                       hl.coalesce(above(ht, ht.xy_het_frac, cohort, 'max_xy_het_frac'), False))
    return clauses


def hail_excess_het(ht, cohort):
    """Hail's excess-het test fails the site: the smallest per-allele p-value is below
    excess_het_p and, if excess_het_f is set, that allele's inbreeding coefficient is below it."""
    f_limit = threshold(ht, cohort, 'excess_het_f')
    clear_excess = hl.if_else(hl.is_defined(f_limit), hl.coalesce(ht.excess_het_f < f_limit, False), True)
    return hl.coalesce(below(ht, ht.p_excess_het, cohort, 'excess_het_p'), False) & clear_excess


def any_of(clauses):
    out = clauses[0]
    for c in clauses[1:]:
        out = out | c
    return out


# AoU: a site passes unless a clause is True (or the result is missing, as in the original)
aou = hl.read_table(args.input_aou_stats_ht_url)
aou_excess = hail_excess_het(aou, 'aou')
if AOU_EXCESS_HET == 'flag':
    aou_excess = hl.if_else(aou.region == 'nonPAR', aou_excess, aou.f_ExcessHet)
aou_fail = any_of(common_clauses(aou, 'aou') + [
    below(aou, aou.max_ac, 'aou', 'min_ac'),
    aou.f_LowQual,
    aou.f_NO_HQ_GENOTYPES,
    aou_excess])
aou.filter(aou_fail, keep=False).select().write(args.output_aou_pass_ht_url, overwrite=True)

# AnVIL: a locus is excluded where a clause is True
anvil = hl.read_table(args.input_anvil_stats_ht_url)
anvil_exclude = any_of(common_clauses(anvil, 'anvil') + [
    hail_excess_het(anvil, 'anvil'),
    anvil.n_hq == 0])
anvil.filter(anvil_exclude).key_by('locus').select().distinct() \
    .write(args.output_anvil_exclusion_ht_url, overwrite=True)

# Kept sites. Both tables are keyed by locus (the pass list by locus and alleles, of
# which locus is a prefix), so this is an ordered join with no shuffle.
pass_ht = hl.read_table(args.output_aou_pass_ht_url).key_by('locus').select()
excl_ht = hl.read_table(args.output_anvil_exclusion_ht_url)
pass_ht.filter(hl.is_missing(excl_ht[pass_ht.locus])).distinct().write(args.output_keep_ht_url, overwrite=True)
keep_ht = hl.read_table(args.output_keep_ht_url)
n_aou, n_pass, n_excl, n_keep = aou.count(), pass_ht.count(), excl_ht.count(), keep_ht.count()
print(f"AoU sites: {n_aou}; pass: {n_pass}; AnVIL excluded loci: {n_excl}; kept: {n_keep}")
if n_keep == 0:
    sys.exit("No site passes the filters; nothing to combine")

intervals = load_intervals(args.bed_url, CHR_VAR)


def trim_padding(rd):
    """Remove what the padding before each interval brought in but nothing needs.

    read_region keeps reference rows that start up to `pad` bp before an interval, for the
    sake of the few blocks that begin there and reach into it. Every other entry in the
    padding is blanked, and rows left with no entry are dropped. Leaving them in is not
    harmless: the combiner partitions by the larger cohort's rows, so a megabase of empty
    AnVIL padding rows lands in one partition and that single task runs for hours.
    """
    ivs = hl.literal(intervals, hl.tarray(hl.tinterval(hl.tlocus('GRCh38'))))
    starts = hl.literal(sorted(iv.start.position for iv in intervals))
    pos = rd.locus.position
    rd = rd.annotate_rows(_inside=ivs.any(lambda iv: iv.contains(rd.locus)),
                          _next_start=hl.min(starts.filter(lambda start: start > pos)))
    # inside an interval: keep. In the padding: keep only blocks that reach the next interval.
    rd = rd.filter_entries(rd._inside | (block_end(rd) >= rd._next_start))
    rd = rd.filter_rows(rd._inside | (hl.agg.count() > 0))
    return rd.drop('_inside', '_next_start')


def filtered_vds(vds_path, label, drop_var, drop_ref, rename_samples):
    rd, vd = read_region(vds_path, CHR_VAR, intervals, args.ref_pad_bp, label)

    # Variant data: kept sites only
    vd = vd.filter_rows(hl.is_defined(keep_ht[vd.locus]))
    vd = drop_present(vd, drop_var)
    vd = vd.select_entries(*SELECT_VAR)

    # Reference data: drop the padding rows that no interval needs (see trim_padding)
    rd = trim_padding(rd)
    rd = drop_present(rd, drop_ref)
    rd = rd.select_entries(*SELECT_REF)

    if rename_samples:
        # str.replace takes a regular expression: '$' matches the end of the ID,
        # so this appends '_C' to every AnVIL sample ID (same call as the original Step 3b).
        vd = vd.key_cols_by(s=vd.s.replace('$', '_C'))
        rd = rd.key_cols_by(s=rd.s.replace('$', '_C'))
    return hl.vds.VariantDataset(rd, vd)


filtered_vds(args.input_aou_vds_url, 'AoU', DROP_AOU_VAR, DROP_AOU_REF, rename_samples=False) \
    .write(args.output_aou_filtered_vds_url, overwrite=True)
filtered_vds(args.input_anvil_vds_url, 'AnVIL', DROP_ANVIL_VAR, DROP_ANVIL_REF, rename_samples=True) \
    .write(args.output_anvil_filtered_vds_url, overwrite=True)

# ---------------------------------------------------------------------- combine
# Same call as the original Step 4. The combiner refuses to run if a finished VDS
# already exists at the output path, so every filter set needs its own filter_tag
# (the WDL puts it in the output names).
combiner = hl.vds.new_combiner(
    output_path=args.output_combined_vds_url,
    temp_path=temp_dir,
    vds_paths=[args.output_aou_filtered_vds_url, args.output_anvil_filtered_vds_url],
    use_genome_default_intervals=True)
combiner.run()

# ----------------------------------------------------------------------- report
end = datetime.now()
time_est = (end - start).total_seconds() / 3600
export_report = (
    f"Chro:\t{CHR_VAR}\n"
    f"Path:\t{args.output_combined_vds_url}\n"
    f"Filt:\t{args.filter_spec}\n"
    f"Time:\t{time_est}\n"
    f"Rigs:\t{args.rig_cost}\n"
    f"Cost:\t{args.rig_cost * time_est}")
local_fname = "tmp_step_report.txt"
with open(local_fname, 'w') as fp:
    fp.write(export_report)
run_subprocess(f"gcloud storage cp {local_fname} {args.output_report_url}",
               f"Error copying report - {local_fname} to {args.output_report_url}")
