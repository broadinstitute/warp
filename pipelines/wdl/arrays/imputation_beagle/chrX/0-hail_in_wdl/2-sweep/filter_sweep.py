#!/usr/bin/env python3
"""Score filter sets offline: precision and recall of SNP genotypes for each candidate
filter set, without re-running Hail, Beagle or vcfdist.

Inputs
  --aou-stats / --anvil-stats   the site-statistics TSVs written by Steps 1 and 2
                                (<contig>.step1.aou.site_stats.tsv.bgz, ...step2.anvil...)
  --vcfdist-dir                 a folder holding the outputs of ONE VcfdistEvaluation run
                                on the unfiltered evaluation callset (the WDL output
                                eval_callset_vcf, run with do_naively_phase = true).
                                Give --gcs-prefix as well to download them there first;
                                an --out-dir already filled by vcfdist_report_sex.py works too.
  --labels                      stratification labels, in vcfdist_bed_files order
  --config NAME:SPEC            a filter set to score (repeatable); SPEC is exactly what
                                the Hail WDL takes as filter_spec. --configs-file holds
                                one NAME<tab>SPEC per line. "original:none" is always scored.

Example
    python3 filter_sweep.py \\
        --aou-stats gs://<out>/chrX.step1.aou.site_stats.tsv.bgz \\
        --anvil-stats gs://<out>/chrX.step2.anvil.site_stats.tsv.bgz \\
        --gcs-prefix gs://<bucket>/submissions/<vcfdist-on-eval-callset> \\
        --vcfdist-dir vcfdist_evalset --labels PAR1,nonPAR \\
        --config par_relaxed:aou.par.min_call_rate=0.8,anvil.par.min_call_rate=0.8 \\
        --config anvil_depth6:anvil.nonpar.min_avsad=6 \\
        --config aou_flag_in_par:aou.excess_het=flag

How a filter set is scored
  A site is kept if its AoU row passes the AoU filters and its position is not on the
  AnVIL exclusion list (the rule Step 3 applies). Every SNP row of vcfdist's query and
  truth tables is assigned to its AoU site. Under a filter set:
      query TP / FP   = TP / FP rows at kept sites
      truth TP        = truth TP rows at kept sites; all other truth rows count as FN
  Precision and recall are computed per sample and averaged, as in evaluation_summary.tsv.

What this does not capture
  * SNPs only (indels cannot be matched to unnormalised VDS rows by position).
  * Genotypes are scored before phasing, and the callset is phased naively, so
    superclusters holding two or more hets can be scored wrong. Use the numbers to rank
    filter sets, then confirm the winner with one full run.
  * vcfdist scores variants in clusters; removing a neighbour can change an outcome.

Needs Python 3.6+, pandas and numpy.
"""
import argparse
import collections
import concurrent.futures
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd

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


# ------------------------------------------------------------------ site tables
def fetch(path, work_dir):
    if not path.startswith("gs://"):
        return path
    os.makedirs(work_dir, exist_ok=True)
    dst = os.path.join(work_dir, path.rstrip("/").rsplit("/", 1)[-1])
    if not os.path.exists(dst):
        print("Downloading %s ..." % path, file=sys.stderr)
        subprocess.check_call(["gcloud", "storage", "cp", path, dst])
    return dst


def read_stats(path, work_dir):
    df = pd.read_csv(fetch(path, work_dir), sep="\t", compression="gzip", dtype=str,
                     keep_default_na=False, na_values=["NA", ""])
    for c in df.columns:
        if c in ("contig", "ref", "alts", "region", "het_test_set") or c.startswith("eh_"):
            continue                                  # eh_*: per-allele lists, not used here
        if c.startswith("f_"):
            # true / false / missing -> 1.0 / 0.0 / NaN
            df[c] = df[c].str.lower().map({"true": 1.0, "false": 0.0})
        else:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df["pos"] = df["pos"].astype(np.int64)
    region = df["region"].fillna("").astype(str)
    df["_nonpar"] = (region == "nonPAR").values
    df["_par"] = region.str.startswith("PAR").values
    return df


def region_thresholds(df, by_region):
    """Per-row threshold array (NaN where the threshold is off)."""
    def val(v):
        return np.nan if v is None else v
    out = np.full(len(df), val(by_region["all"]), dtype=float)
    out[df["_nonpar"].values] = val(by_region["nonpar"])
    out[df["_par"].values] = val(by_region["par"])
    return out


class Clauses:
    """Three-valued OR of filter clauses, as Hail evaluates them: each clause is True,
    False or missing per row."""
    def __init__(self, n):
        self.true = np.zeros(n, dtype=bool)
        self.missing = np.zeros(n, dtype=bool)

    def add(self, true, missing=None):
        self.true |= true
        if missing is not None:
            self.missing |= missing


def compare(x, thr, op):
    """(true, missing) for x < thr or x > thr; an off threshold (NaN) gives False, not missing."""
    on = ~np.isnan(thr)
    x_missing = np.isnan(x)
    with np.errstate(invalid="ignore"):
        hit = (x < thr) if op == "<" else (x > thr)
    return on & ~x_missing & hit, on & x_missing


def common_clauses(df, thresholds, cl):
    not_dense = np.isnan(df["call_rate"].values)        # site skipped by the dense pass
    t, _ = compare(df["mean_sum_ad"].values, region_thresholds(df, thresholds["min_avsad"]), "<")
    cl.add((df["n_ad"].values >= 1) & t)
    t, m = compare(df["call_rate"].values, region_thresholds(df, thresholds["min_call_rate"]), "<")
    cl.add(t, m)
    t, m = compare(df["mean_gq"].values, region_thresholds(df, thresholds["min_gq"]), "<")
    cl.add(t, m & not_dense)                            # NaN mean GQ on a densified site is not "missing"
    if "xy_het_frac" in df.columns:
        t, _ = compare(df["xy_het_frac"].values, region_thresholds(df, thresholds["max_xy_het_frac"]), ">")
        # only where at least xy_het_min_carriers XY samples carry a non-reference allele (same rule as Step 3)
        min_carriers = region_thresholds(df, thresholds["xy_het_min_carriers"])
        with np.errstate(invalid="ignore"):
            enough = np.isnan(min_carriers) | (df["n_xy_nonref"].values >= min_carriers)
        cl.add(df["_nonpar"].values & enough & t)


def hail_excess_het(df, thresholds):
    """Smallest per-allele p below excess_het_p and, if excess_het_f is set, that allele's
    inbreeding coefficient below it (same rule as Step 3)."""
    t, _ = compare(df["p_excess_het"].values, region_thresholds(df, thresholds["excess_het_p"]), "<")
    f_limit = region_thresholds(df, thresholds["excess_het_f"])
    f_hit, _ = compare(df["excess_het_f"].values, f_limit, "<")
    return t & (np.isnan(f_limit) | f_hit)


# ------------------------------------------------- excess-het tests of the earlier runs
def excess_het_p_exact(n_hom_ref, n_het, n_hom_var, alpha):
    """True where Hail's one-sided exact Hardy-Weinberg test (excess heterozygotes) gives
    p < alpha: p = P(at least n_het heterozygotes | sample size, allele counts)
    (Levene-Haldane distribution, as in hl.hardy_weinberg_test(..., one_sided=True))."""
    from scipy.special import gammaln
    n_hom_ref, n_het, n_hom_var = (np.nan_to_num(np.asarray(x, dtype=float)) for x in (n_hom_ref, n_het, n_hom_var))
    n = n_hom_ref + n_het + n_hom_var
    n_a = n_het + 2 * np.minimum(n_hom_ref, n_hom_var)          # minor allele count
    n_b = 2 * n - n_a
    out = np.zeros(len(n), dtype=bool)
    with np.errstate(invalid="ignore", divide="ignore"):
        expected = n_a * n_b / (2 * n - 1)
        q = n_a / (2 * n)
        sd = np.sqrt(4 * n * (q * (1 - q)) ** 2)                # never larger than the true one
        z = (n_het - expected) / sd
    # far from significance (p >= ~1e-2 at z = 2): decided without the exact sum
    for i in np.flatnonzero((n > 0) & (n_het > expected) & ((z > 2) | (alpha > 1e-3))):
        k = np.arange(n_het[i] % 2, n_a[i] + 1, 2)              # possible heterozygote counts
        logp = (gammaln(n[i] + 1) - gammaln((n_a[i] - k) / 2 + 1) - gammaln(k + 1) - gammaln((n_b[i] - k) / 2 + 1)
                + k * np.log(2.0) + gammaln(n_a[i] + 1) + gammaln(n_b[i] + 1) - gammaln(2 * n[i] + 1))
        prob = np.exp(logp - logp.max())
        out[i] = prob[k >= n_het[i]].sum() / prob.sum() < alpha
    return out


def first_number(col):
    return pd.to_numeric(col.astype(str).str.split(",").str[0], errors="coerce").values


def all_sample_excess_het(df, alpha, xy_ref_counted=True):
    """The excess-het filter of the earlier runs, rebuilt from the per-site counts:
    hl.variant_qc(...).p_value_excess_het < alpha on the densified data with every sample in.
      - biallelic sites only (variant_qc gives no p-value at multiallelic sites);
      - a haploid non-reference call counts as a heterozygote (AC[1] - 2 * hom_var);
      - xy_ref_counted: whether XY reference calls in non-PAR count as hom-ref (they do when
        diploid 0/0; a haploid 0 call is not a homozygote for variant_qc).
    Samples without a sex label are not in the non-PAR counts (well under 1% of samples)."""
    het, hom = first_number(df["eh_n_het"]), first_number(df["eh_n_hom_var"])
    ref = df["n_test_diploid"].values - het - hom
    xx_only = df["_nonpar"].values & (df["het_test_set"].astype(str).str.lower() != "all").values
    xy = {c: np.nan_to_num(df[c].values) for c in ("n_xy_called", "n_xy_nonref", "n_xy_het", "n_xy_haploid_nonref")}
    het = het + np.where(xx_only, xy["n_xy_het"] + xy["n_xy_haploid_nonref"], 0)
    hom = hom + np.where(xx_only, xy["n_xy_nonref"] - xy["n_xy_het"] - xy["n_xy_haploid_nonref"], 0)
    if xy_ref_counted:
        ref = ref + np.where(xx_only, xy["n_xy_called"] - xy["n_xy_nonref"], 0)
    flagged = excess_het_p_exact(np.maximum(ref, 0), het, np.maximum(hom, 0), alpha)
    return flagged & (df["n_alleles"].values == 2)


def check_exact_test(df, name, alpha=3.4e-6):
    """Where the table's own test used every sample, the rebuilt test must agree with it."""
    use = (df["n_alleles"].values == 2) & ~(df["_nonpar"].values) & ~np.isnan(df["p_excess_het"].values)
    if not use.any():
        return
    d = df[use]
    het, hom = first_number(d["eh_n_het"]), first_number(d["eh_n_hom_var"])
    mine = excess_het_p_exact(d["n_test_diploid"].values - het - hom, het, hom, alpha)
    table = d["p_excess_het"].values < alpha
    print("check of the rebuilt excess-het test, %s biallelic sites outside non-PAR: %d sites, %d flagged in the "
          "table, %d by the rebuilt test, %d disagree" % (name, use.sum(), table.sum(), mine.sum(), (mine != table).sum()))


EARLIER_RUNS = collections.OrderedDict([
    # name: (filter spec for the thresholds, where AoU uses its precomputed flag, XY refs counted in AoU non-PAR)
    ("pre_handoff_run", ("none", "everywhere", None)),
    ("post_report_run_xyref_as_homref", ("both.par.min_call_rate=0.8,both.par.min_avsad=10,both.par.min_gq=25", "par", True)),
    ("post_report_run_xyref_not_counted", ("both.par.min_call_rate=0.8,both.par.min_avsad=10,both.par.min_gq=25", "par", False)),
])


def earlier_run_site_kept(aou, anvil, name):
    """Sites kept by the filters of an earlier run.
    pre_handoff_run:  original thresholds; AoU: precomputed ExcessHet flag everywhere;
                      AnVIL: all-sample variant_qc test.
    post_report_run:  PAR thresholds relaxed (depth 10, call rate 0.8, GQ 25); AoU: precomputed flag in
                      PAR, all-sample variant_qc test in non-PAR; AnVIL: all-sample variant_qc test."""
    spec, flag_where, xy_ref = EARLIER_RUNS[name]
    thresholds, _ = parse_filter_spec(spec)
    alpha = FILTER_DEFAULTS["aou"]["excess_het_p"]
    flag = aou["f_ExcessHet"].values
    if flag_where == "everywhere":
        aou_true, aou_missing = flag == 1.0, np.isnan(flag)
    else:
        nonpar = aou["_nonpar"].values
        aou_true = np.where(nonpar, all_sample_excess_het(aou, alpha, xy_ref), flag == 1.0)
        aou_missing = ~nonpar & np.isnan(flag)
    passes = aou_pass(aou, thresholds["aou"], "hail", excess=(aou_true, aou_missing))
    excluded = anvil_excluded_positions(anvil, thresholds["anvil"], excess=all_sample_excess_het(anvil, alpha))
    return passes, excluded


def aou_pass(aou, thresholds, mode, excess=None):
    """Boolean array: AoU rows that pass. A row is removed if any clause is True, or if
    none is True and one is missing (Hail's filter drops rows whose predicate is missing).
    excess: (true, missing) arrays that replace the excess-het clause (earlier runs)."""
    cl = Clauses(len(aou))
    common_clauses(aou, thresholds, cl)
    t, m = compare(aou["max_ac"].values, region_thresholds(aou, thresholds["min_ac"]), "<")
    cl.add(t, m)
    for flag in ("f_LowQual", "f_NO_HQ_GENOTYPES"):
        v = aou[flag].values
        cl.add(v == 1.0, np.isnan(v))
    if excess is not None:
        cl.add(*excess)
        return ~(cl.true | cl.missing)
    excess = hail_excess_het(aou, thresholds)
    if mode == "flag":
        v = aou["f_ExcessHet"].values
        outside = ~aou["_nonpar"].values
        cl.add(np.where(outside, v == 1.0, excess), outside & np.isnan(v))
    else:
        cl.add(excess)
    return ~(cl.true | cl.missing)


def anvil_excluded_positions(anvil, thresholds, excess=None):
    """Positions of AnVIL rows where any clause is True."""
    cl = Clauses(len(anvil))
    common_clauses(anvil, thresholds, cl)
    cl.add(hail_excess_het(anvil, thresholds) if excess is None else excess)
    cl.add(anvil["n_hq"].values == 0)
    return np.unique(anvil["pos"].values[cl.true])


def trim(pos, ref, alt):
    """Minimal representation of one allele: drop shared trailing, then leading, bases."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def snp_index(df):
    """{(pos, ref, alt): row number} for every SNP allele of every row, after trimming."""
    index = {}
    for i, (pos, ref, alts) in enumerate(zip(df["pos"].tolist(), df["ref"].tolist(), df["alts"].tolist())):
        if len(ref) == 1 and len(alts) == 1:
            index.setdefault((pos, ref, alts), i)
            continue
        for alt in str(alts).split(","):
            p, r, a = trim(pos, ref, alt)
            if len(r) == 1 and len(a) == 1:
                index.setdefault((p, r, a), i)
    return index


# -------------------------------------------------------------- vcfdist outputs
KINDS = {".query.tsv": "query", ".truth.tsv": "truth"}
SINGLETONS = ("sample_manifest.tsv",)


def run(cmd):
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if r.returncode != 0:
        sys.exit("command failed: %s\n%s" % (" ".join(cmd), r.stderr.strip()))
    return r.stdout


def download(prefix, out_dir, billing_project):
    extra = ["--billing-project", billing_project] if billing_project else []
    listing = run(["gcloud", "storage", "ls", *extra, prefix.rstrip("/") + "/**"])
    paths = [p.strip() for p in listing.splitlines() if p.startswith("gs://")]
    paths = [p for p in paths if p.rsplit("/", 1)[-1] in SINGLETONS or any(p.endswith(s) for s in KINDS)]
    jobs = []
    for p in paths:
        local = os.path.join(out_dir, "raw", p[len(prefix.rstrip("/")) + 1:])
        if not os.path.exists(local):
            os.makedirs(os.path.dirname(local), exist_ok=True)
            jobs.append((p, local))
    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        list(pool.map(lambda j: run(["gcloud", "storage", "cp", *extra, j[0], j[1]]), jobs))
    print("%d vcfdist files found, %d downloaded" % (len(paths), len(jobs)), file=sys.stderr)


def stratum_index(rel_path):
    """Index of the outer scatter (the stratification) from a Cromwell path."""
    before, _, after = rel_path.partition("call-Vcfdist")
    outer = re.findall(r"shard-(\d+)", before)
    if outer:
        return int(outer[-1])
    inner = re.findall(r"shard-(\d+)", after)
    return int(inner[0]) if len(inner) >= 2 else 0


def attempt_rank(rel_path):
    m = re.search(r"attempt-(\d+)", rel_path)
    return int(m.group(1)) if m else 1


def index_files(out_dir):
    """{(stratum_idx, sample, kind): path}, keeping the latest attempt."""
    files, singles = {}, {}
    root = os.path.join(out_dir, "raw")
    for d, _, names in os.walk(root):
        for name in names:
            path = os.path.join(d, name)
            rel = os.path.relpath(path, root)
            if name in SINGLETONS:
                singles[name] = path
                continue
            for sfx, kind in KINDS.items():
                if name.endswith(sfx):
                    key = (stratum_index(rel), name[:-len(sfx)], kind)
                    if key not in files or attempt_rank(rel) >= attempt_rank(os.path.relpath(files[key], root)):
                        files[key] = path
    return files, singles


def load_variants(files, label, sex):
    """One row per SNP row of every query / truth table: stratum, sample, sex, side, pos, ref, alt, err."""
    frames = []
    for (j, sample, side), path in sorted(files.items()):
        v = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        need = {"POS", "REF", "ALT", "TYPE", "ERRTYPE"}
        if not need <= set(v.columns):
            sys.exit("%s lacks columns %s" % (path, sorted(need - set(v.columns))))
        v = v[v["TYPE"] == "SNP"]
        frames.append(pd.DataFrame({
            "stratum": label.get(j, "stratum%d" % j), "sample": sample, "sex": sex.get(sample, "NA"), "side": side,
            "pos": pd.to_numeric(v["POS"], errors="coerce").fillna(-1).astype(np.int64).values,
            "ref": v["REF"].values, "alt": v["ALT"].values, "err": v["ERRTYPE"].values}))
    return pd.concat(frames, ignore_index=True)


# ----------------------------------------------------------------------- score
def score(variants, kept):
    """Mean per-sample SNP counts, precision and recall for one filter set, by stratum and sex."""
    v = variants.assign(kept=kept)
    q, t = v[v["side"] == "query"], v[v["side"] == "truth"]
    per = pd.DataFrame({
        "QUERY_TP": q[(q["err"] == "TP") & q["kept"]].groupby(["stratum", "sex", "sample"]).size(),
        "QUERY_FP": q[(q["err"] == "FP") & q["kept"]].groupby(["stratum", "sex", "sample"]).size(),
        "TRUTH_TP": t[(t["err"] == "TP") & t["kept"]].groupby(["stratum", "sex", "sample"]).size(),
        "TRUTH_ALL": t[t["err"].isin(["TP", "FN"])].groupby(["stratum", "sex", "sample"]).size(),
    }).fillna(0.0).reset_index()
    per["TRUTH_FN"] = per["TRUTH_ALL"] - per["TRUTH_TP"]
    per["PREC"] = per["QUERY_TP"] / (per["QUERY_TP"] + per["QUERY_FP"]).where(lambda x: x > 0)
    per["RECALL"] = per["TRUTH_TP"] / per["TRUTH_ALL"].where(lambda x: x > 0)
    per["F1"] = 2 * per["PREC"] * per["RECALL"] / (per["PREC"] + per["RECALL"]).where(lambda x: x > 0)
    both = pd.concat([per, per.assign(sex="all")], ignore_index=True)
    cols = ["QUERY_TP", "QUERY_FP", "TRUTH_TP", "TRUTH_FN", "PREC", "RECALL", "F1"]
    out = both.groupby(["stratum", "sex"])[cols].mean()
    out.insert(0, "samples", both.groupby(["stratum", "sex"])["sample"].nunique())
    return out.reset_index()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aou-stats", required=True)
    ap.add_argument("--anvil-stats", required=True)
    ap.add_argument("--vcfdist-dir", required=True, help="local folder for the vcfdist outputs (files under <dir>/raw)")
    ap.add_argument("--gcs-prefix", help="gs:// path of the VcfdistEvaluation submission to download into --vcfdist-dir")
    ap.add_argument("--billing-project")
    ap.add_argument("--labels", default="", help="comma-separated stratification labels, in vcfdist_bed_files order")
    ap.add_argument("--sex-tsv", help="local two-column sample/sex file; default is the run's sample_manifest.tsv")
    ap.add_argument("--config", action="append", default=[], help="NAME:SPEC (repeatable)")
    ap.add_argument("--configs-file", help="one NAME<tab>SPEC per line")
    ap.add_argument("--earlier-runs", action="store_true",
                    help="also score the filters of the pre-handoff run and the post-report run, rebuilt from "
                         "the site statistics (rows pre_handoff_run and post_report_run_*)")
    ap.add_argument("--out", default="filter_sweep.tsv")
    args = ap.parse_args()

    configs = collections.OrderedDict([("original", "none")])
    for item in args.config:
        if ":" not in item:
            ap.error("--config must be NAME:SPEC, got %r" % item)
        name, spec = item.split(":", 1)
        configs[name] = spec
    if args.configs_file:
        for line in open(args.configs_file):
            f = line.split()
            if len(f) >= 2 and not line.startswith("#"):
                configs[f[0]] = f[1]
    parsed = collections.OrderedDict((name, parse_filter_spec(spec)) for name, spec in configs.items())   # fail early

    if args.gcs_prefix:
        download(args.gcs_prefix, args.vcfdist_dir, args.billing_project)
    files, singles = index_files(args.vcfdist_dir)
    if not files:
        sys.exit("no vcfdist query/truth tables under %s/raw" % args.vcfdist_dir)
    labels = [x for x in args.labels.split(",") if x]
    label = {j: (labels[j] if j < len(labels) else "stratum%d" % j) for j in sorted({k[0] for k in files})}
    sex = {}
    if args.sex_tsv:
        sex = dict(line.split()[:2] for line in open(args.sex_tsv) if len(line.split()) >= 2)
    elif "sample_manifest.tsv" in singles:
        m = pd.read_csv(singles["sample_manifest.tsv"], sep="\t", dtype=str)
        if {"eval_sample", "sex"} <= set(m.columns):
            sex = dict(zip(m["eval_sample"], m["sex"]))

    work = os.path.join(args.vcfdist_dir, "site_stats")
    print("Reading site statistics ...", file=sys.stderr)
    aou, anvil = read_stats(args.aou_stats, work), read_stats(args.anvil_stats, work)
    index = snp_index(aou)
    aou_pos = aou["pos"].values
    print("AoU sites: %d   AnVIL sites: %d" % (len(aou), len(anvil)))

    variants = load_variants(files, label, sex)
    other = variants[~variants["err"].isin(["TP", "FP", "FN"])]
    if len(other):
        print("ignoring %d rows with other outcomes: %s" % (len(other), dict(other["err"].value_counts())))

    # vcfdist positions may be 0-based: pick the offset under which matched query SNPs are found
    tp = variants[(variants["side"] == "query") & (variants["err"] == "TP")].drop_duplicates(["pos", "ref", "alt"]).head(20000)
    rates = {off: sum((p + off, r, a) in index for p, r, a in zip(tp["pos"], tp["ref"], tp["alt"])) / max(1, len(tp))
             for off in (0, 1)}
    offset = max(rates, key=rates.get)
    print("query TP SNPs found in the AoU table with position offset 0: %.3f, +1: %.3f -> using %+d"
          % (rates[0], rates[1], offset))
    row = np.array([index.get((p + offset, r, a), -1)
                    for p, r, a in zip(variants["pos"], variants["ref"], variants["alt"])], dtype=np.int64)
    in_table = row >= 0
    for side in ("query", "truth"):
        s = variants["side"].values == side
        print("%s SNP rows: %d, of which %d (%.1f%%) are not in the AoU site table and can never be kept"
              % (side, s.sum(), (s & ~in_table).sum(), 100.0 * (s & ~in_table).sum() / max(1, s.sum())))
    if rates[offset] < 0.9:
        print("WARNING: fewer than 90% of matched query SNPs were found in the AoU table; "
              "check that the statistics and the evaluation callset come from the same region and run")

    results = []
    ceiling = score(variants, in_table)
    ceiling.insert(0, "config", "(no filters)")
    ceiling.insert(1, "sites_kept", len(aou))
    results.append(ceiling)
    for name, (thresholds, mode) in parsed.items():
        passes = aou_pass(aou, thresholds["aou"], mode)
        excluded = anvil_excluded_positions(anvil, thresholds["anvil"])
        site_kept = passes & ~np.isin(aou_pos, excluded)
        kept = in_table & site_kept[np.where(in_table, row, 0)]
        res = score(variants, kept)
        res.insert(0, "config", name)
        res.insert(1, "sites_kept", int(site_kept.sum()))
        results.append(res)
        print("%-24s AoU pass %d, AnVIL excluded positions %d, kept sites %d   [%s]"
              % (name, passes.sum(), len(excluded), site_kept.sum(), configs[name]))

    if args.earlier_runs:
        check_exact_test(aou, "AoU")
        check_exact_test(anvil, "AnVIL")
        for name in EARLIER_RUNS:
            passes, excluded = earlier_run_site_kept(aou, anvil, name)
            site_kept = passes & ~np.isin(aou_pos, excluded)
            kept = in_table & site_kept[np.where(in_table, row, 0)]
            res = score(variants, kept)
            res.insert(0, "config", name)
            res.insert(1, "sites_kept", int(site_kept.sum()))
            results.append(res)
            print("%-34s AoU pass %d, AnVIL excluded positions %d, kept sites %d"
                  % (name, passes.sum(), len(excluded), site_kept.sum()))

    table = pd.concat(results, ignore_index=True)
    table.to_csv(args.out, sep="\t", index=False, float_format="%.4f")
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 2000)
    for (stratum, sx), d in table.groupby(["stratum", "sex"], sort=True):
        print("\n== %s, sex %s (mean per sample, SNPs) ==" % (stratum, sx))
        print(d.drop(columns=["stratum", "sex"]).round(4).to_string(index=False))
    print("\nWrote %s" % args.out)


if __name__ == "__main__":
    main()
