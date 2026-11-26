#!/usr/bin/env python3
"""Download the outputs of a VcfdistEvaluation submission and summarise them.

Example (Terra terminal or notebook):

    python3 vcfdist_report.py \\
        --gcs-prefix gs://<bucket>/submissions/<submission-id> \\
        --labels PAR1,nonPAR

--labels must be in the same order as the workflow's vcfdist_bed_files input.
Needs Python 3.6+ and pandas (matplotlib optional). Re-running skips files that are already downloaded. To re-analyse without
touching the bucket, pass --out-dir alone (omit --gcs-prefix).

Outputs, all in --out-dir:
    report.txt              everything below as text: paste this back
    per_sample.tsv          one row per stratum and sample
    variant_breakdown.tsv   counts by stratum, sex, side, type, zygosity, outcome
    truth_snp_lookup.tsv    with --site-tables: every truth SNP and why it is or is not in the panel
    position_bins.tsv       counts per position bin, summed over samples
    *.png                   plots of the same tables
"""
import argparse
import collections
import concurrent.futures
import os
import re
import subprocess
import sys
import time

import numpy as np
import pandas as pd

KINDS = {
    ".precision-recall-summary.tsv": "pr_summary",
    ".phasing-summary.tsv": "phasing_summary",
    ".switchflips.tsv": "switchflips",
    ".phase-blocks.tsv": "phase_blocks",
    ".superclusters.tsv": "superclusters",
    ".query.tsv": "query",
    ".truth.tsv": "truth",
}
SINGLETONS = ("evaluation_summary.tsv", "sample_manifest.tsv")


def log(msg):
    print(msg, file=sys.stderr, flush=True)


# ----------------------------------------------------------------- download

def run(cmd):
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if r.returncode != 0:
        sys.exit("command failed: %s\n%s" % (" ".join(cmd), r.stderr.strip()))
    return r.stdout


def list_bucket(prefix, billing_project):
    extra = ["--billing-project", billing_project] if billing_project else []
    out = run(["gcloud", "storage", "ls", *extra, prefix.rstrip("/") + "/**"])
    return [line.strip() for line in out.splitlines() if line.startswith("gs://")]


def wanted(path):
    name = path.rsplit("/", 1)[-1]
    return name in SINGLETONS or any(name.endswith(sfx) for sfx in KINDS)


def download(prefix, out_dir, billing_project):
    log("Listing %s ..." % prefix)
    paths = [p for p in list_bucket(prefix, billing_project) if wanted(p)]
    print("%d result files found under %s" % (len(paths), prefix), flush=True)
    extra = ["--billing-project", billing_project] if billing_project else []
    jobs = []
    for p in paths:
        local = os.path.join(out_dir, "raw", p[len(prefix.rstrip("/")) + 1:])
        if not os.path.exists(local):
            os.makedirs(os.path.dirname(local), exist_ok=True)
            jobs.append((p, local))
    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        list(pool.map(lambda j: run(["gcloud", "storage", "cp", *extra, j[0], j[1]]), jobs))
    print("%d downloaded, %d already present" % (len(jobs), len(paths) - len(jobs)), flush=True)


# ------------------------------------------------------------------- index

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
    if not os.path.isdir(root):
        return files, singles
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


# ----------------------------------------------------------------- helpers

def read(path):
    try:
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except Exception as exc:
        print("could not read %s: %s" % (path, exc), file=sys.stderr)
        return pd.DataFrame()


def col(df, *names):
    """First column of df whose name matches one of names (case-insensitive)."""
    lower = {c.lower(): c for c in df.columns}
    for n in names:
        if n.lower() in lower:
            return lower[n.lower()]
    return None


def num(series):
    return pd.to_numeric(series, errors="coerce")


def show(df, **kw):
    return "(empty)" if df is None or len(df) == 0 else df.to_string(**kw)


def with_all_sex(df, sex_col="sex"):
    """Append a copy of df with sex='all' so XX, XY, and all can be grouped together."""
    all_df = df.copy()
    all_df[sex_col] = "all"
    return pd.concat([df, all_df], ignore_index=True)


# ----------------------------------------------------- missed-site lookup

SITE_TABLES = {"aou": "aou.sites.tsv.bgz", "anvil": "anvil.sites.tsv.bgz",
               "step1": "step1_pass.sites.tsv.bgz", "excl": "anvil_exclusion.sites.tsv.bgz"}


def trim(pos, ref, alt):
    """Minimal representation of one allele: drop shared trailing, then leading, bases."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def load_site_tables(location, out_dir, billing_project):
    """Read the tables written by export_site_tables.py from a local or gs:// folder."""
    local = location
    if location.startswith("gs://"):
        local = os.path.join(out_dir, "site_tables")
        os.makedirs(local, exist_ok=True)
        extra = ["--billing-project", billing_project] if billing_project else []
        for name in SITE_TABLES.values():
            dst = os.path.join(local, name)
            if not os.path.exists(dst):
                log("Downloading %s ..." % name)
                r = subprocess.run(["gcloud", "storage", "cp", *extra, location.rstrip("/") + "/" + name, dst],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                if r.returncode != 0:
                    if os.path.exists(dst):
                        os.remove(dst)
                    log("site table %s not available" % name)
    tables = {}
    for key, name in SITE_TABLES.items():
        path = os.path.join(local, name)
        if os.path.exists(path):
            log("Reading %s ..." % path)
            tables[key] = pd.read_csv(path, sep="\t", compression="gzip", dtype=str, keep_default_na=False)
            log("  loaded %s (%d rows)" % (name, len(tables[key])))
    return tables


def snp_index(df):
    """{(pos, ref, alt): row number} for every SNP allele of every row, after trimming."""
    index = {}
    positions = pd.to_numeric(df["pos"].astype(str).str.strip("[]\"'"), errors="coerce").fillna(0).astype(int).tolist()
    refs = df["ref"].astype(str).str.strip("[]\"'").tolist()
    alts_col = df["alts"].astype(str).str.strip("[]\"'").tolist()

    for i, (pos, ref, alts) in enumerate(zip(positions, refs, alts_col)):
        if len(ref) == 1 and len(alts) == 1:
            index.setdefault((pos, ref, alts), i)
            continue
        for alt in alts.split(","):
            p, r, a = trim(pos, ref, alt)
            if len(r) == 1 and len(a) == 1:
                index.setdefault((p, r, a), i)
    return index


def precompute_cohort_arrays(df, fallback_n, fallback_xx_frac):
    """Precompute vectorized per-row XX/XY/all counts and AF arrays for fast dict lookup."""
    pos = pd.to_numeric(df["pos"].astype(str).str.strip("[]\"'"), errors="coerce").fillna(0).astype(int).values
    is_par = (pos <= 2781479) | (pos >= 155701383)

    def _num(col_name, default_val=0.0):
        if col_name in df.columns:
            return pd.to_numeric(df[col_name].astype(str).str.strip("[]\"'"), errors="coerce").fillna(default_val).values
        return np.full(len(df), default_val, dtype=np.float64)

    n_all = _num("n_samples_all", float(fallback_n))
    n_xx = _num("n_samples_xx", float(int(round(fallback_n * fallback_xx_frac))))
    n_xy = _num("n_samples_xy", float(fallback_n - int(round(fallback_n * fallback_xx_frac))))

    ac_all = _num("ac_all") if "ac_all" in df.columns else _num("cohort_max_AC")
    het_xx = _num("n_dip_het_xx") if "n_dip_het_xx" in df.columns else _num("n_dip_het")
    hom_xx = _num("n_dip_hom_var_xx") if "n_dip_hom_var_xx" in df.columns else _num("n_dip_hom_var")
    ac_xx = _num("ac_xx") if "ac_xx" in df.columns else (het_xx + 2.0 * hom_xx)

    hap_xy = _num("n_hap_nonref_xy") if "n_hap_nonref_xy" in df.columns else _num("n_hap_nonref")
    ac_xy = _num("ac_xy") if "ac_xy" in df.columns else hap_xy

    has_hap = (_num("n_hap_nonref") > 0).any()
    an_all = np.where(is_par, 2.0 * n_all, np.where(has_hap, 2.0 * n_xx + n_xy, 2.0 * n_all))
    an_xx = 2.0 * n_xx
    an_xy = np.where(is_par, 2.0 * n_xy, np.where(has_hap, n_xy, 2.0 * n_xy))
    carriers_xx = het_xx + hom_xx

    return {
        "ac_all": ac_all.astype(np.int64),
        "ac_xx": ac_xx.astype(np.int64),
        "ac_xy": ac_xy.astype(np.int64),
        "af_all": np.where(an_all > 0, ac_all / an_all, np.nan),
        "af_xx": np.where(an_xx > 0, ac_xx / an_xx, np.nan),
        "af_xy": np.where(an_xy > 0, ac_xy / an_xy, np.nan),
        "xx_carrier_het_share": np.where(carriers_xx > 0, het_xx / carriers_xx, np.nan),
    }


class SiteLookup:
    def __init__(self, tables):
        aou = tables["aou"]
        log("Indexing AoU sites...")
        self.aou_index = snp_index(aou)
        self.aou_pos = pd.to_numeric(aou["pos"], errors="coerce").fillna(0).astype(int).values
        flag = lambda c: (aou[c].str.lower() == "true").values if c in aou.columns else np.zeros(len(aou), dtype=bool)
        self.excess_het = flag("f_ExcessHet")
        self.other_flag = flag("f_LowQual") | flag("f_NO_HQ_GENOTYPES")
        self.max_ac = pd.to_numeric(aou["cohort_max_AC"], errors="coerce").fillna(0).values
        self.aou_stats = precompute_cohort_arrays(aou, 414830, 0.606)

        if "anvil" in tables:
            log("Indexing AnVIL sites...")
            self.anvil_index = snp_index(tables["anvil"])
            self.anvil_stats = precompute_cohort_arrays(tables["anvil"], 100749, 0.455)
        else:
            self.anvil_index = None
            self.anvil_stats = None

        self.step1_index = snp_index(tables["step1"]) if "step1" in tables else None
        self.excl_pos = set(pd.to_numeric(tables["excl"]["pos"], errors="coerce").dropna().astype(int)) if "excl" in tables else None
        self.offset = 0

    def category(self, pos, ref, alt):
        key = (pos + self.offset, ref, alt)
        i = self.aou_index.get(key)
        if i is None:
            if self.anvil_index is not None and key in self.anvil_index:
                return "1 not in AoU VDS (in AnVIL only)"
            return "1 not in AoU VDS"
        # A site in the Step 1 output passed Step 1 whatever flags it carries
        # (Step 1 no longer filters on the AoU ExcessHet flag in non-PAR).
        passed_step1 = self.step1_index is not None and key in self.step1_index
        if not passed_step1:
            if self.other_flag[i]:
                return "3 AoU LowQual / NO_HQ flag"
            if self.max_ac[i] < 2:
                return "4 AoU allele count < 2"
            if self.excess_het[i]:
                return "2 AoU ExcessHet flag (may not be the clause that removed it)"
            if self.step1_index is None:
                return "5 in AoU, unflagged (give the Step 1 table to resolve)"
            return "5 AoU other Step 1 clause (call rate / depth / GQ / Hail excess het)"
        if self.excl_pos is not None and int(self.aou_pos[i]) in self.excl_pos:
            return "6 passed Step 1, on AnVIL exclusion list"
        return "7 site kept by both steps"

    def enrich_truth_df(self, ts):
        """Attach sex-stratified AoU and AnVIL counts/AFs to truth_snps via O(1) dict lookups."""
        n = len(ts)
        keys = [(int(p) + self.offset, r, a) for p, r, a in zip(ts["pos"], ts["ref"], ts["alt"])]

        aou_rows = [self.aou_index.get(k, -1) for k in keys]
        aou_idx = np.asarray(aou_rows, dtype=np.int64)
        has_aou = aou_idx >= 0
        ts["aou_present"] = has_aou
        for col_name, arr in self.aou_stats.items():
            out_col = np.full(n, np.nan if arr.dtype.kind == "f" else 0, dtype=arr.dtype)
            out_col[has_aou] = arr[aou_idx[has_aou]]
            ts["aou_" + col_name] = np.round(out_col, 6) if arr.dtype.kind == "f" else out_col

        if self.anvil_index is not None and self.anvil_stats is not None:
            anvil_rows = [self.anvil_index.get(k, -1) for k in keys]
            anvil_idx = np.asarray(anvil_rows, dtype=np.int64)
            has_anvil = anvil_idx >= 0
            ts["anvil_present"] = has_anvil
            for col_name, arr in self.anvil_stats.items():
                out_col = np.full(n, np.nan if arr.dtype.kind == "f" else 0, dtype=arr.dtype)
                out_col[has_anvil] = arr[anvil_idx[has_anvil]]
                ts["anvil_" + col_name] = np.round(out_col, 6) if arr.dtype.kind == "f" else out_col
        return ts


# ------------------------------------------------------------------ report

def main():
    t0 = time.time()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gcs-prefix", help="gs:// path of the submission (or workflow) directory")
    ap.add_argument("--out-dir", default="vcfdist_results")
    ap.add_argument("--labels", default="", help="comma-separated stratification labels, in vcfdist_bed_files order")
    ap.add_argument("--sex-tsv", help="local two-column sample/sex file; default is the run's sample_manifest.tsv")
    ap.add_argument("--bin-size", type=int, default=100000)
    ap.add_argument("--billing-project")
    ap.add_argument("--site-tables", help="folder (local or gs://) written by export_site_tables.py; "
                                          "adds section 7, which says why each missed truth SNP is absent")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    if args.gcs_prefix:
        download(args.gcs_prefix, args.out_dir, args.billing_project)

    log("Indexing vcfdist result files under %s/raw ..." % args.out_dir)
    files, singles = index_files(args.out_dir)
    if not files:
        sys.exit("no vcfdist result files under %s/raw" % args.out_dir)
    log("Found %d files across %d strata and %d samples (%.1fs)" % (
        len(files), len({k[0] for k in files}), len({k[1] for k in files}), time.time() - t0
    ))

    labels = [x for x in args.labels.split(",") if x]
    strata = sorted({k[0] for k in files})
    label = {j: (labels[j] if j < len(labels) else "stratum%d" % j) for j in strata}
    samples = sorted({k[1] for k in files})

    sex = {}
    if args.sex_tsv:
        for line in open(args.sex_tsv):
            f = line.split()
            if len(f) >= 2:
                sex[f[0]] = f[1]
    elif "sample_manifest.tsv" in singles:
        m = read(singles["sample_manifest.tsv"])
        if col(m, "eval_sample") and col(m, "sex"):
            sex = dict(zip(m[col(m, "eval_sample")], m[col(m, "sex")]))

    out = []
    def section(title):
        out.append("\n" + "=" * 78 + "\n" + title + "\n" + "=" * 78)

    # ---- 0. inventory
    section("0. FILES")
    inv = collections.Counter((label[k[0]], k[2]) for k in files)
    out.append("samples: %d   strata: %s" % (len(samples), ", ".join("%d=%s" % (j, label[j]) for j in strata)))
    out.append(show(pd.Series(inv).unstack(fill_value=0)))
    out.append("sex known for %d of %d samples: %s" % (sum(s in sex for s in samples), len(samples),
                                                      dict(collections.Counter(sex.get(s, "NA") for s in samples))))

    # ---- 1. schemas
    section("1. SCHEMAS (first stratum, first sample that has each file)")
    for kind in sorted(set(KINDS.values())):
        keys = sorted(k for k in files if k[2] == kind)
        if not keys:
            continue
        df = read(files[keys[0]])
        out.append("\n-- %s  (%d rows)  columns: %s" % (kind, len(df), list(df.columns)))
        out.append(show(df.head(4), index=False))
    if "evaluation_summary.tsv" in singles:
        out.append("\n-- evaluation_summary.tsv")
        out.append(open(singles["evaluation_summary.tsv"]).read().strip())

    # ---- per-sample pass
    log("Parsing per-sample vcfdist outputs...")
    rows, variants, bins, block_rows, phase_values = [], [], [], [], collections.Counter()
    truth_snps = []
    for j in strata:
        for s in samples:
            row = {"stratum": label[j], "sample": s, "sex": sex.get(s, "NA")}

            if (j, s, "pr_summary") in files:
                pr = read(files[(j, s, "pr_summary")])
                vt, thr = col(pr, "VAR_TYPE", "TYPE"), col(pr, "THRESHOLD", "MIN_QUAL")
                if vt and thr:
                    for _, r in pr[pr[thr] == "NONE"].iterrows():
                        for c in ("TRUTH_TP", "QUERY_TP", "TRUTH_FN", "QUERY_FP", "PREC", "RECALL"):
                            if c in pr.columns:
                                row["%s_%s" % (r[vt], c)] = pd.to_numeric(r[c], errors="coerce")

            if (j, s, "phasing_summary") in files:
                ps = read(files[(j, s, "phasing_summary")])
                if len(ps):
                    for c in ps.columns:
                        row["PHASING_" + c] = pd.to_numeric(ps[c].iloc[0], errors="coerce")

            if (j, s, "switchflips") in files:
                sf = read(files[(j, s, "switchflips")])
                tcol, pcol = col(sf, "SWITCH_TYPE", "TYPE", "ERR_TYPE"), col(sf, "START", "POS")
                row["switchflip_rows"] = len(sf)
                if tcol:
                    for value, n in sf[tcol].value_counts().items():
                        row["switchflips_" + str(value)] = n
                if pcol and len(sf):
                    d = pd.DataFrame({"pos": num(sf[pcol]), "what": sf[tcol] if tcol else "ERROR"}).dropna()
                    d["bin"] = (d["pos"] // args.bin_size * args.bin_size).astype(int)
                    for (bin_start, value), n in d.groupby(["bin", "what"]).size().items():
                        bins.append((label[j], row["sex"], bin_start, "switchflips_" + str(value), n))

            if (j, s, "phase_blocks") in files:
                pb = read(files[(j, s, "phase_blocks")])
                row["phase_block_rows"] = len(pb)
                a, b = col(pb, "START"), col(pb, "STOP", "END")
                size = num(pb[col(pb, "SIZE")]) if col(pb, "SIZE") else (num(pb[b]) - num(pb[a]) if a and b else None)
                if size is not None and len(size):
                    row["phase_block_median_bp"] = size.median()
                    row["phase_block_max_bp"] = size.max()
                    for v in size.dropna():
                        block_rows.append((label[j], s, v))

            if (j, s, "superclusters") in files:
                sc = read(files[(j, s, "superclusters")])
                row["supercluster_rows"] = len(sc)
                for pcol in (col(sc, "SC_PHASE", "PHASE"), col(sc, "PHASE_STATE")):
                    if pcol:
                        for value, n in sc[pcol].value_counts().items():
                            row["%s=%s" % (pcol, value)] = n
                            phase_values[(label[j], pcol, str(value))] += n
                pcol = col(sc, "SC_PHASE", "PHASE")
                if pcol:
                    row["sc_informative"] = int((~sc[pcol].isin(["?", "-", ".", ""])).sum())

            for side in ("truth", "query"):
                if (j, s, side) not in files:
                    continue
                v = read(files[(j, s, side)])
                pos, vt = col(v, "POS", "START"), col(v, "TYPE", "VAR_TYPE")
                et, hap = col(v, "ERRTYPE", "ERR_TYPE", "ERROR_TYPE", "BD"), col(v, "HAP", "HAPLOTYPE")
                if not (pos and et):
                    continue
                site = [c for c in (col(v, "CONTIG", "CHROM"), pos, col(v, "REF"), col(v, "ALT")) if c]
                g = v.groupby(site, sort=False)
                per_site = pd.DataFrame({
                    "pos": num(g[pos].first()),
                    "type": g[vt].first() if vt else "NA",
                    "outcome": g[et].agg(lambda x: "/".join(sorted(set(x)))),
                    "zyg": g[hap].nunique().map({1: "het", 2: "hom"}).fillna("other") if hap else "NA",
                })
                row["%s_sites" % side] = len(per_site)
                if side == "truth" and args.site_tables and len(site) == 4:
                    for key, r in per_site[per_site["type"] == "SNP"].iterrows():
                        truth_snps.append((label[j], row["sex"], s, int(r["pos"]), key[2], key[3], r["zyg"], r["outcome"]))
                row["%s_het_TP_sites" % side] = int(((per_site["zyg"] == "het") & (per_site["outcome"] == "TP")).sum())
                row["%s_pos_min" % side], row["%s_pos_max" % side] = per_site["pos"].min(), per_site["pos"].max()
                for (t, z, o), n in per_site.groupby(["type", "zyg", "outcome"]).size().items():
                    variants.append((label[j], row["sex"], s, side, t, z, o, n))
                binned = per_site.dropna(subset=["pos"]).copy()
                binned["bin"] = (binned["pos"] // args.bin_size * args.bin_size).astype(int)
                for (bin_start, o), n in binned.groupby(["bin", "outcome"]).size().items():
                    bins.append((label[j], row["sex"], bin_start, "%s_%s" % (side, o), n))
            rows.append(row)
    log("Parsed per-sample outputs (%.1fs)" % (time.time() - t0))

    per_sample = pd.DataFrame(rows)
    sw, fl, nb = "PHASING_SWITCH_ERRORS", "PHASING_FLIP_ERRORS", "PHASING_PHASE_BLOCKS"
    if {sw, fl, nb} <= set(per_sample.columns):
        for unit, name in (("sc_informative", "per_informative_supercluster"), ("query_het_TP_sites", "per_het_TP_site")):
            if unit in per_sample.columns:
                opp = (per_sample[unit] - per_sample[nb]).where(lambda x: x > 0)
                per_sample["switch_rate_" + name] = per_sample[sw] / opp
                per_sample["switchflip_rate_" + name] = (per_sample[sw] + per_sample[fl]) / opp
    per_sample.to_csv(os.path.join(args.out_dir, "per_sample.tsv"), sep="\t", index=False)
    numeric = [c for c in per_sample.columns if c not in ("stratum", "sample", "sex")]

    section("2. PER SAMPLE")
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 200)
    out.append(show(per_sample.round(4), index=False))

    section("3. MEAN PER SAMPLE, BY STRATUM AND SEX (XX, XY, all)")
    numeric = [c for c in numeric if pd.api.types.is_numeric_dtype(per_sample[c])]
    ps_with_all = with_all_sex(per_sample)
    out.append(show(ps_with_all.groupby(["stratum", "sex"])[numeric].mean().round(3).T))

    section("4. VARIANT SITES BY TYPE, ZYGOSITY AND OUTCOME (mean count per sample)")
    vb = pd.DataFrame(variants, columns=["stratum", "sex", "sample", "side", "type", "zygosity", "outcome", "n"])
    vb.to_csv(os.path.join(args.out_dir, "variant_breakdown.tsv"), sep="\t", index=False)
    if len(vb):
        vb_with_all = with_all_sex(vb)
        n_samples = ps_with_all.groupby(["stratum", "sex"])["sample"].nunique()
        tot = vb_with_all.groupby(["stratum", "sex", "side", "type", "zygosity", "outcome"])["n"].sum()
        mean = (tot / [n_samples[(k[0], k[1])] for k in tot.index]).round(1)
        out.append(show(mean.unstack("outcome", fill_value=0)))
    else:
        out.append("(no per-variant columns recognised; see the query/truth schemas in section 1)")

    section("5. COUNTS PER %d bp BIN, SUMMED OVER SAMPLES" % args.bin_size)
    pbins = pd.DataFrame(bins, columns=["stratum", "sex", "bin_start", "what", "n"])
    wide = None
    if len(pbins):
        wide = pbins.groupby(["stratum", "bin_start", "what"])["n"].sum().unstack("what", fill_value=0)
        wide.to_csv(os.path.join(args.out_dir, "position_bins.tsv"), sep="\t")
        out.append(show(wide))
        out.append("\nSame, by sex (XX, XY, all — totals over all bins):")
        pb_with_all = with_all_sex(pbins)
        out.append(show(pb_with_all.groupby(["stratum", "sex", "what"])["n"].sum().unstack("what", fill_value=0)))

    section("6. PHASING")
    keep = [c for c in per_sample.columns if c.startswith(("PHASING_PHASE", "PHASING_SWITCH_ERR", "PHASING_FLIP", "phase_block_rows",
                                                           "supercluster", "SC_PHASE", "PHASE_STATE", "sc_informative",
                                                           "query_het_TP", "truth_het_TP", "switch_rate", "switchflip_rate"))]
    out.append(show(per_sample[["stratum", "sample", "sex"] + keep].round(5), index=False))
    rate_cols = [c for c in per_sample.columns if c.startswith(("switch_rate", "switchflip_rate"))]
    if rate_cols:
        out.append("\nPooled over samples (sum of errors / sum of opportunities), by stratum and sex (XX, XY, all):")
        pooled = []
        for (st, sx), d in ps_with_all.groupby(["stratum", "sex"]):
            rec = {"stratum": st, "sex": sx, "samples": len(d), "switch": d[sw].sum(), "flip": d[fl].sum()}
            for unit in ("sc_informative", "query_het_TP_sites"):
                if unit in d.columns:
                    opp = (d[unit] - d[nb]).clip(lower=0).sum()
                    rec[unit + "_opportunities"] = opp
                    rec["switch_rate/" + unit] = d[sw].sum() / opp if opp else float("nan")
                    rec["switchflip_rate/" + unit] = (d[sw].sum() + d[fl].sum()) / opp if opp else float("nan")
            pooled.append(rec)
        out.append(show(pd.DataFrame(pooled).round(5), index=False))
    out.append("\nSupercluster phase values, summed over samples: %s" % dict(phase_values))
    blocks = pd.DataFrame(block_rows, columns=["stratum", "sample", "size_bp"])
    if len(blocks):
        out.append("\nPhase block sizes (bp), all samples:")
        out.append(show(blocks.groupby("stratum")["size_bp"].describe().round(0)))

    if args.site_tables:
        log("Running Section 7 site-table lookup...")
        section("7. WHY TRUTH SNPs ARE (NOT) IN THE PANEL")
        tables = load_site_tables(args.site_tables, args.out_dir, args.billing_project)
        if "aou" not in tables:
            out.append("aou.sites.tsv.bgz not found in %s" % args.site_tables)
        elif not truth_snps:
            out.append("no truth SNPs collected (the truth.tsv columns were not recognised)")
        else:
            lookup = SiteLookup(tables)
            ts = pd.DataFrame(truth_snps, columns=["stratum", "sex", "sample", "pos", "ref", "alt", "zygosity", "outcome"])
            tp = ts[ts["outcome"] == "TP"].drop_duplicates(["pos", "ref", "alt"]).head(20000)
            rates = {}
            for off in (0, 1):
                rates[off] = sum((p + off, r, a) in lookup.aou_index for p, r, a in zip(tp["pos"], tp["ref"], tp["alt"])) / max(1, len(tp))
            lookup.offset = max(rates, key=rates.get)
            out.append("site tables found: %s" % ", ".join(sorted(tables)))
            out.append("matched (TP) truth SNPs found in the AoU table with position offset 0: %.3f, offset +1: %.3f -> using %+d"
                       % (rates[0], rates[1], lookup.offset))
            ts["category"] = [lookup.category(p, r, a) for p, r, a in zip(ts["pos"], ts["ref"], ts["alt"])]
            ts = lookup.enrich_truth_df(ts)
            ts.to_csv(os.path.join(args.out_dir, "truth_snp_lookup.tsv"), sep="\t", index=False)

            ts["missed"] = ts["outcome"].map(lambda o: "matched (TP)" if o == "TP" else "missed (%s)" % o)
            ts_with_all = with_all_sex(ts)
            n_s = ts_with_all.groupby(["stratum", "sex"])["sample"].nunique()
            out.append("\nMean truth SNP sites per sample, by reason (XX, XY, all):")
            cnt = ts_with_all.groupby(["stratum", "sex", "missed", "category"]).size()
            out.append(show((cnt / [n_s[(k[0], k[1])] for k in cnt.index]).round(1).unstack(["stratum", "sex", "missed"], fill_value=0)))

            out.append("\nShare of MISSED truth SNPs by reason (XX, XY, all):")
            miss = ts_with_all[ts_with_all["outcome"] != "TP"]
            share = miss.groupby(["stratum", "sex", "category"]).size().unstack(["stratum", "sex"], fill_value=0)
            out.append(show((share / share.sum()).round(3)))

            out.append("\nSame, split by zygosity in the truth:")
            miss_orig = ts[ts["outcome"] != "TP"]
            share_zyg = miss_orig.groupby(["stratum", "zygosity", "category"]).size().unstack(["stratum", "zygosity"], fill_value=0)
            out.append(show((share_zyg / share_zyg.sum()).round(3)))

            aou = tables["aou"]
            if "f_ExcessHet" in aou.columns and "n_dip_het" in aou.columns:
                a = aou[["pos", "f_ExcessHet"]].copy()
                for c in ("n_dip_het", "n_dip_hom_var", "n_hap_nonref", "cohort_max_AC",
                          "n_dip_het_xx", "n_dip_hom_var_xx", "n_dip_het_xy", "n_hap_nonref_xy"):
                    if c in aou.columns:
                        a[c] = pd.to_numeric(aou[c], errors="coerce")
                a["pos"] = a["pos"].astype(int)
                a = a[a["cohort_max_AC"] >= 1000]
                a["region"] = ["PAR1" if p <= 2781479 else "nonPAR" for p in a["pos"]]
                het_xx = a["n_dip_het_xx"] if "n_dip_het_xx" in a.columns else a["n_dip_het"]
                hom_xx = a["n_dip_hom_var_xx"] if "n_dip_hom_var_xx" in a.columns else a["n_dip_hom_var"]
                hap_xy = a["n_hap_nonref_xy"] if "n_hap_nonref_xy" in a.columns else a["n_hap_nonref"]
                dip_het_xy = a["n_dip_het_xy"] if "n_dip_het_xy" in a.columns else 0

                a["xx_carrier_het_share"] = het_xx / (het_xx + hom_xx).where(lambda x: x > 0)
                a["xy_haploid_share"] = hap_xy / (hap_xy + dip_het_xy).where(lambda x: x > 0)
                a["AC_bin"] = pd.cut(a["cohort_max_AC"], [1e3, 1e4, 5e4, 1e5, 2e5, 4e5, 1e7], right=False)
                g = a.groupby(["region", "AC_bin", "f_ExcessHet"])
                out.append("\nAoU sites with cohort allele count >= 1000: medians by ExcessHet flag")
                out.append(show(pd.DataFrame({"sites": g.size(),
                                              "xx_carrier_het_share": g["xx_carrier_het_share"].median().round(3),
                                              "xy_haploid_share": g["xy_haploid_share"].median().round(3)})
                                .query("sites > 0")))

    report = "\n".join(out)
    with open(os.path.join(args.out_dir, "report.txt"), "w") as fh:
        fh.write(report + "\n")
    print(report, flush=True)

    # ---- plots
    log("Generating summary plots...")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib not available: skipping plots", file=sys.stderr)
        return

    colors = {v: c for v, c in zip(sorted(per_sample["sex"].unique()), ["#1f77b4", "#d62728", "#7f7f7f", "#2ca02c"])}
    if {"SNP_PREC", "SNP_RECALL"} <= set(per_sample.columns):
        fig, axes = plt.subplots(1, len(strata), figsize=(5 * len(strata), 4.5), squeeze=False)
        for ax, j in zip(axes[0], strata):
            d = per_sample[per_sample["stratum"] == label[j]]
            for sx, dd in d.groupby("sex"):
                ax.scatter(dd["SNP_RECALL"], dd["SNP_PREC"], label=sx, color=colors[sx], s=45)
            ax.set(title="%s: SNP precision vs recall per sample" % label[j], xlabel="recall", ylabel="precision",
                   xlim=(0, 1.02), ylim=(0, 1.02))
            ax.grid(alpha=0.3)
            ax.legend(title="sex")
        fig.tight_layout()
        fig.savefig(os.path.join(args.out_dir, "precision_recall_per_sample.png"), dpi=130)
        plt.close(fig)

    if wide is not None:
        for j in strata:
            if label[j] not in wide.index.get_level_values(0):
                continue
            d = wide.loc[label[j]]
            d = d.loc[:, (d != 0).any()]
            fig, axes = plt.subplots(len(d.columns), 1, figsize=(11, 1.7 * len(d.columns) + 1), sharex=True, squeeze=False)
            for ax, c in zip(axes[:, 0], d.columns):
                ax.bar(d.index.astype(float) / 1e6, d[c], width=args.bin_size / 1e6 * 0.9, align="edge", color="#444444")
                ax.set_ylabel(c, rotation=0, ha="right", va="center", fontsize=8)
                ax.grid(alpha=0.3)
            axes[-1, 0].set_xlabel("position (Mb)")
            fig.suptitle("%s: counts per %d kb bin, summed over samples" % (label[j], args.bin_size // 1000))
            fig.tight_layout()
            fig.savefig(os.path.join(args.out_dir, "position_bins_%s.png" % label[j]), dpi=130)
            plt.close(fig)

    if len(blocks):
        fig, ax = plt.subplots(figsize=(7, 4))
        for st, d in blocks.groupby("stratum"):
            ax.hist(d["size_bp"].clip(lower=1), bins=40, alpha=0.6, label=st)
        ax.set(xlabel="phase block size (bp)", ylabel="blocks", title="Phase block sizes, all samples")
        ax.set_xscale("log")
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(args.out_dir, "phase_block_sizes.png"), dpi=130)
        plt.close(fig)
    print("\nWrote report.txt, TSVs and plots to %s (%.1fs)" % (args.out_dir, time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
