#!/usr/bin/env python3
"""Plot full-cohort Hardy-Weinberg deviation, H_obs/H_exp ratio, and signed
-log10(p_HWE) stratified by XX, XY, and all across all sites evaluated in
Steps 1 and 2, highlighting:
  1) lrWGS truth membership,
  2) AoU / Step2 ExcessHet flag status, and
  3) Secondary Step1 / Step2 QC filter status (f_other_qc).

Usage:
    python3 plot_hw_full_cohort.py \
        --aou-sites gs://<bucket>/site_tables/aou.sites.tsv.bgz \
        --anvil-sites gs://<bucket>/site_tables/anvil.sites.tsv.bgz \
        --truth-lookup vcfdist_results/truth_snp_lookup.tsv \
        --out-dir hw_plots_full_cohort
"""

import argparse
import math
import os
import shutil
import subprocess
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PVAL_CLIP = 50.0
EH_CUTOFF_P = 3.4e-6
EH_CUTOFF_NLP = -math.log10(EH_CUTOFF_P)
LN10 = math.log(10.0)


def trim(pos, ref, alt):
    """Minimal representation of one allele: drop shared trailing, then leading, bases."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def parse_numeric_col(series):
    """Parse plain numbers or Hail-exported bracketed strings like '[123]'."""
    cleaned = series.fillna("0").astype(str).str.strip().str.strip("[]\"'")
    has_comma = cleaned.str.contains(",", regex=False)
    if has_comma.any():
        def _max_val(val):
            nums = [float(x) for x in val.split(",") if x.strip() not in ("", "NA", "nan", "None")]
            return max(nums) if nums else 0.0
        out = pd.to_numeric(cleaned.where(~has_comma, "0"), errors="coerce").fillna(0.0)
        out.loc[has_comma] = cleaned.loc[has_comma].map(_max_val)
        return out
    return pd.to_numeric(cleaned, errors="coerce").fillna(0.0)


def parse_bool_col(df, col_name):
    """Parse Hail boolean column ('true'/'false') with False default if missing."""
    if col_name not in df.columns:
        return pd.Series(False, index=df.index, dtype=bool)
    return df[col_name].fillna("").astype(str).str.strip().str.lower() == "true"


def _exact_hwe_scalar(n_aa, n_ab, n_bb, cache):
    """Memoized exact two-sided HWE -log10(p) (Wigginton et al. 2005) for small counts."""
    n = n_aa + n_ab + n_bb
    if n <= 0:
        return np.nan
    n_a = 2 * n_aa + n_ab
    n_b = 2 * n_bb + n_ab
    n_rare = min(n_a, n_b)
    if n_rare == 0:
        return 0.0

    key = (n, n_rare, n_ab)
    if key in cache:
        return cache[key]

    n_common = 2 * n - n_rare
    mid = int(round(n_rare * n_common / (2.0 * n)))
    if (mid % 2) != (n_rare % 2):
        mid = mid + 1 if mid < n_rare else mid - 1
    probs = {mid: 1.0}

    curr_hets = mid
    curr_homr = (n_rare - mid) // 2
    curr_homc = (n_common - mid) // 2
    val = 1.0
    while curr_hets <= n_rare - 2:
        val *= (4.0 * curr_homr * curr_homc) / ((curr_hets + 2.0) * (curr_hets + 1.0))
        curr_hets += 2
        curr_homr -= 1
        curr_homc -= 1
        probs[curr_hets] = val

    curr_hets = mid
    curr_homr = (n_rare - mid) // 2
    curr_homc = (n_common - mid) // 2
    val = probs[mid]
    while curr_hets >= 2:
        val *= (curr_hets * (curr_hets - 1.0)) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
        curr_hets -= 2
        curr_homr += 1
        curr_homc += 1
        probs[curr_hets] = val

    total = sum(probs.values())
    p_obs = probs.get(n_ab, 0.0) / total
    pval = sum(v for v in probs.values() if (v / total) <= p_obs * (1.0 + 1e-12)) / total
    res = -math.log10(max(pval, 1e-300))
    cache[key] = res
    return res


def signed_hwe_neg_log10_p_vec(n_aa, n_ab, n_bb):
    """Vectorized signed -log10(p_HWE) across NumPy arrays."""
    n_aa = np.asarray(n_aa, dtype=np.float64)
    n_ab = np.asarray(n_ab, dtype=np.float64)
    n_bb = np.asarray(n_bb, dtype=np.float64)

    n = n_aa + n_ab + n_bb
    n_a = 2.0 * n_aa + n_ab
    n_b = 2.0 * n_bb + n_ab
    n_rare = np.minimum(n_a, n_b)

    out = np.full(len(n), np.nan, dtype=np.float64)
    valid = np.isfinite(n) & (n > 0) & (n_rare > 0)
    if not np.any(valid):
        return out

    large = valid & (n_rare >= 400)
    if np.any(large):
        nl = n[large]
        p = n_a[large] / (2.0 * nl)
        q = 1.0 - p
        e_aa = nl * p * p
        e_ab = 2.0 * nl * p * q
        e_bb = nl * q * q
        x2 = (
            ((n_aa[large] - e_aa) ** 2) / e_aa
            + ((n_ab[large] - e_ab) ** 2) / e_ab
            + ((n_bb[large] - e_bb) ** 2) / e_bb
        )
        z = np.sqrt(np.maximum(0.0, x2) * 0.5)
        t = 1.0 / (1.0 + 0.3275911 * z)
        poly = t * (
            0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + 1.061405429 * t)))
        )
        nlp = np.maximum(0.0, (z * z - np.log(np.maximum(poly, 1e-300))) / LN10)
        sign = np.where(n_ab[large] >= e_ab, 1.0, -1.0)
        out[large] = sign * nlp

    small_idx = np.where(valid & (n_rare < 400))[0]
    if len(small_idx):
        cache = {}
        for idx in small_idx:
            aa = int(round(n_aa[idx]))
            ab = int(round(n_ab[idx]))
            bb = int(round(n_bb[idx]))
            nd = aa + ab + bb
            p = (ab + 2.0 * bb) / (2.0 * nd)
            e_ab = nd * 2.0 * p * (1.0 - p)
            nlp = _exact_hwe_scalar(aa, ab, bb, cache)
            out[idx] = (1.0 if ab >= e_ab else -1.0) * nlp

    return out


def fetch_local(path, work_dir, local_name):
    if path.startswith("gs://"):
        dst = os.path.join(work_dir, local_name)
        subprocess.run(["gcloud", "storage", "cp", path, dst], check=True)
        return dst
    return path


def load_site_table(path, work_dir, local_name, fallback_n, fallback_xx_frac):
    """Load one cohort's per-allele site table (including multiallelic sites and Step1/2 QC flags)."""
    local_path = fetch_local(path, work_dir, local_name)
    raw_df = pd.read_csv(local_path, sep="\t", compression="gzip", dtype=str)
    print(f"  {local_name}: {len(raw_df):,} raw allele rows", flush=True)

    df = raw_df.copy()
    num_cols = [
        c for c in df.columns
        if c.startswith(("pos", "n_alleles", "n_samples_", "site_max_AC", "cohort_max_AC", "ac_", "n_dip_", "n_hap_"))
    ]
    for col in num_cols:
        df[col] = parse_numeric_col(df[col])

    df["ref"] = df["ref"].fillna("").astype(str).str.strip().str.strip("[]\"'")
    df["alts"] = df["alts"].fillna("").astype(str).str.strip().str.strip("[]\"'")

    site_ac_col = "site_max_AC" if "site_max_AC" in df.columns else "cohort_max_AC"
    df = df[(df["alts"] != "") & (df[site_ac_col] >= 2) & (df["cohort_max_AC"] >= 1)].copy()

    if df["alts"].str.contains(",", regex=False).any():
        print(f"  WARNING: {local_name} contains unexploded multiallelic rows; filtering to biallelic.", flush=True)
        df = df[~df["alts"].str.contains(",", regex=False)].copy()

    if df.empty:
        raise ValueError(f"{local_name} had 0 valid rows!")

    trimmed = [trim(int(p), r, a) for p, r, a in zip(df["pos"], df["ref"], df["alts"])]
    df["pos"] = [t[0] for t in trimmed]
    df["ref"] = [t[1] for t in trimmed]
    df["alt"] = [t[2] for t in trimmed]
    df = df.drop_duplicates(["pos", "ref", "alt"]).copy()

    df["region"] = np.where(
        df["pos"] <= 2781479, "PAR1",
        np.where(df["pos"] >= 155701383, "PAR2", "nonPAR"),
    )

    req_sex_cols = {
        "n_samples_xx", "n_samples_xy",
        "n_dip_het_xx", "n_dip_hom_var_xx",
        "n_dip_het_xy", "n_dip_hom_var_xy", "n_hap_nonref_xy",
    }
    has_sex_cols = req_sex_cols.issubset(df.columns)

    if has_sex_cols:
        n_xx = int(df["n_samples_xx"].max())
        n_xy = int(df["n_samples_xy"].max())
        nonpar_xy_is_haploid = (df.loc[df["region"] == "nonPAR", "n_hap_nonref_xy"] > 0).any()

        xy_nonpar_alt_chroms = np.where(
            nonpar_xy_is_haploid,
            df["n_hap_nonref_xy"],
            df["n_dip_hom_var_xy"] + 0.5 * df["n_dip_het_xy"],
        )

        out = pd.DataFrame({
            "pos": df["pos"],
            "ref": df["ref"],
            "alt": df["alt"],
            "region": df["region"],
            "has_sex_cols": True,
            "nonpar_xy_is_haploid": bool(nonpar_xy_is_haploid),
            "n_xx": n_xx,
            "n_xy": n_xy,
            "n_dip_het_xx": df["n_dip_het_xx"],
            "n_dip_hom_var_xx": df["n_dip_hom_var_xx"],
            "n_dip_het_xy": df["n_dip_het_xy"],
            "n_dip_hom_var_xy": df["n_dip_hom_var_xy"],
            "xy_nonpar_alt_chroms": xy_nonpar_alt_chroms,
            "n_dip_het_all": df["n_dip_het_xx"] + df["n_dip_het_xy"],
            "n_dip_hom_var_all": df["n_dip_hom_var_xx"] + df["n_dip_hom_var_xy"],
        })
    else:
        n_xx = int(round(fallback_n * fallback_xx_frac))
        n_xy = int(fallback_n) - n_xx
        nonpar_xy_is_haploid = (df.loc[df["region"] == "nonPAR", "n_hap_nonref"] > 0).any()
        xy_nonpar_alt_chroms = np.where(nonpar_xy_is_haploid, df["n_hap_nonref"], 0.0)

        out = pd.DataFrame({
            "pos": df["pos"],
            "ref": df["ref"],
            "alt": df["alt"],
            "region": df["region"],
            "has_sex_cols": False,
            "nonpar_xy_is_haploid": bool(nonpar_xy_is_haploid),
            "n_xx": n_xx,
            "n_xy": n_xy,
            "n_dip_het_xx": np.nan,
            "n_dip_hom_var_xx": np.nan,
            "n_dip_het_xy": np.nan,
            "n_dip_hom_var_xy": np.nan,
            "xy_nonpar_alt_chroms": xy_nonpar_alt_chroms,
            "n_dip_het_all": df["n_dip_het"],
            "n_dip_hom_var_all": df["n_dip_hom_var"],
        })

    if "f_ExcessHet" in df.columns:
        out["excess_het"] = parse_bool_col(df, "f_ExcessHet")
    else:
        out["excess_het"] = np.nan

    out["f_call_rate"] = parse_bool_col(df, "f_call_rate")
    out["f_gq"] = parse_bool_col(df, "f_gq")
    out["f_avsad"] = parse_bool_col(df, "f_avsad")
    out["f_outside_bed"] = parse_bool_col(df, "f_outside_bed")
    if "f_lowqual_nohq" in df.columns:
        out["f_lowqual_nohq"] = parse_bool_col(df, "f_lowqual_nohq")
    else:
        out["f_lowqual_nohq"] = parse_bool_col(df, "f_LowQual") | parse_bool_col(df, "f_NO_HQ_GENOTYPES")

    if "f_other_qc" in df.columns:
        out["f_other_qc"] = parse_bool_col(df, "f_other_qc")
    else:
        out["f_other_qc"] = (
            out["f_call_rate"] | out["f_gq"] | out["f_avsad"] | out["f_lowqual_nohq"] | out["f_outside_bed"]
        )
    return out


def score_site_table(df):
    """Compute AF, HW deviation, sex-aware H_obs/H_exp ratio, and signed -log10(p)
    for XX, XY, and all."""
    out = df.copy()
    is_par = out["region"].isin(["PAR1", "PAR2"])
    has_sex = bool(out["has_sex_cols"].iloc[0])
    xy_dip_nonpar = ~bool(out["nonpar_xy_is_haploid"].iloc[0])

    # --- XX (always diploid) ---
    out["n_dip_hom_ref_xx"] = np.maximum(0.0, out["n_xx"] - out["n_dip_het_xx"] - out["n_dip_hom_var_xx"])
    ac_xx = out["n_dip_het_xx"] + 2.0 * out["n_dip_hom_var_xx"]
    an_xx = 2.0 * out["n_xx"]
    out["af_xx"] = np.where(an_xx > 0, ac_xx / an_xx, np.nan)
    out["het_frac_xx"] = np.where(out["n_xx"] > 0, out["n_dip_het_xx"] / out["n_xx"], np.nan)
    exp_xx = 2.0 * out["af_xx"] * (1.0 - out["af_xx"])
    out["hw_dev_xx"] = out["het_frac_xx"] - exp_xx
    out["het_ratio_xx"] = np.where(exp_xx > 0, out["het_frac_xx"] / exp_xx, np.nan)

    # --- XY (diploid in PAR; 1 X chromosome per male in nonPAR) ---
    out["n_dip_hom_ref_xy"] = np.where(
        is_par,
        np.maximum(0.0, out["n_xy"] - out["n_dip_het_xy"] - out["n_dip_hom_var_xy"]),
        np.nan,
    )
    ac_xy = np.where(
        is_par,
        out["n_dip_het_xy"] + 2.0 * out["n_dip_hom_var_xy"],
        out["xy_nonpar_alt_chroms"],
    )
    an_xy = np.where(is_par, 2.0 * out["n_xy"], 1.0 * out["n_xy"])
    out["af_xy"] = np.where(an_xy > 0, ac_xy / an_xy, np.nan)
    out["het_frac_xy"] = np.where(out["n_xy"] > 0, out["n_dip_het_xy"] / out["n_xy"], np.nan)
    exp_xy = 2.0 * out["af_xy"] * (1.0 - out["af_xy"])
    out["hw_dev_xy"] = out["het_frac_xy"] - exp_xy
    out["het_ratio_xy"] = np.where(exp_xy > 0, out["het_frac_xy"] / exp_xy, np.nan)

    # --- all (combined XX + XY) ---
    n_all = out["n_xx"] + out["n_xy"]
    if has_sex:
        ac_all = ac_xx + ac_xy
        an_all = an_xx + an_xy
    else:
        ac_all = np.where(
            is_par | xy_dip_nonpar,
            out["n_dip_het_all"] + 2.0 * out["n_dip_hom_var_all"],
            out["n_dip_het_all"] + 2.0 * out["n_dip_hom_var_all"] + out["xy_nonpar_alt_chroms"],
        )
        an_all = np.where(is_par | xy_dip_nonpar, 2.0 * n_all, 2.0 * out["n_xx"] + out["n_xy"])

    out["af_all"] = np.where(an_all > 0, ac_all / an_all, np.nan)
    out["het_frac_all"] = np.where(n_all > 0, out["n_dip_het_all"] / n_all, np.nan)
    exp_all_2pq = 2.0 * out["af_all"] * (1.0 - out["af_all"])
    out["hw_dev_all"] = out["het_frac_all"] - exp_all_2pq

    dip_share_all = np.where(is_par, 1.0, out["n_xx"] / n_all)
    exp_all_sex_aware = dip_share_all * exp_all_2pq
    out["het_ratio_all"] = np.where(exp_all_sex_aware > 0, out["het_frac_all"] / exp_all_sex_aware, np.nan)

    # --- Vectorized exact HWE p-values on diploid genotypes ---
    out["signed_p_xx"] = signed_hwe_neg_log10_p_vec(
        out["n_dip_hom_ref_xx"], out["n_dip_het_xx"], out["n_dip_hom_var_xx"]
    )
    p_xy_par = signed_hwe_neg_log10_p_vec(
        out["n_dip_hom_ref_xy"], out["n_dip_het_xy"], out["n_dip_hom_var_xy"]
    )
    out["signed_p_xy"] = np.where(
        is_par,
        p_xy_par,
        np.where(out["af_xy"].notna() & (out["af_xy"] > 0), 0.0, np.nan),
    )

    n_hom_ref_all_par = np.maximum(0.0, n_all - out["n_dip_het_all"] - out["n_dip_hom_var_all"])
    p_all_par = signed_hwe_neg_log10_p_vec(
        n_hom_ref_all_par, out["n_dip_het_all"], out["n_dip_hom_var_all"]
    )
    if has_sex:
        out["signed_p_all"] = np.where(is_par, p_all_par, out["signed_p_xx"])
    else:
        p_all_legacy = signed_hwe_neg_log10_p_vec(
            np.maximum(0.0, out["n_xx"] - out["n_dip_het_all"] - out["n_dip_hom_var_all"]),
            out["n_dip_het_all"],
            out["n_dip_hom_var_all"],
        )
        out["signed_p_all"] = np.where(is_par, p_all_par, p_all_legacy)

    if out["excess_het"].isna().any():
        ratio_check = np.where(is_par | (~has_sex), out["het_ratio_all"] > 1.0, out["het_ratio_xx"] > 1.0)
        calc_eh = ratio_check & (out["signed_p_all"] > EH_CUTOFF_NLP)
        out["excess_het"] = out["excess_het"].fillna(calc_eh).astype(bool)
    else:
        out["excess_het"] = out["excess_het"].astype(bool)

    return out


def combine_aou_anvil(aou_df, anvil_df):
    """Merge AoU and AnVIL per-allele counts and Step1/2 flags into AoU+AnVIL."""
    both_have_sex = bool(aou_df["has_sex_cols"].iloc[0]) and bool(anvil_df["has_sex_cols"].iloc[0])
    sum_cols = [
        "n_dip_het_xx", "n_dip_hom_var_xx",
        "n_dip_het_xy", "n_dip_hom_var_xy", "xy_nonpar_alt_chroms",
        "n_dip_het_all", "n_dip_hom_var_all",
    ]
    flag_cols = ["excess_het", "f_call_rate", "f_gq", "f_avsad", "f_lowqual_nohq", "f_outside_bed", "f_other_qc"]
    m = pd.merge(
        aou_df[["pos", "ref", "alt", "region"] + flag_cols + sum_cols],
        anvil_df[["pos", "ref", "alt", "region"] + flag_cols + sum_cols],
        on=["pos", "ref", "alt", "region"],
        how="outer",
        suffixes=("_aou", "_anvil"),
    )

    if both_have_sex:
        for c in sum_cols:
            m[c] = m[f"{c}_aou"].fillna(0.0) + m[f"{c}_anvil"].fillna(0.0)
    else:
        for c in ("n_dip_het_xx", "n_dip_hom_var_xx", "n_dip_het_xy", "n_dip_hom_var_xy"):
            m[c] = np.nan
        for c in ("xy_nonpar_alt_chroms", "n_dip_het_all", "n_dip_hom_var_all"):
            m[c] = m[f"{c}_aou"].fillna(0.0) + m[f"{c}_anvil"].fillna(0.0)

    m["n_xx"] = int(aou_df["n_xx"].iloc[0]) + int(anvil_df["n_xx"].iloc[0])
    m["n_xy"] = int(aou_df["n_xy"].iloc[0]) + int(anvil_df["n_xy"].iloc[0])
    m["has_sex_cols"] = both_have_sex
    m["nonpar_xy_is_haploid"] = True
    for fc in flag_cols:
        m[fc] = m[f"{fc}_aou"].fillna(False).astype(bool) | m[f"{fc}_anvil"].fillna(False).astype(bool)

    keep_cols = ["pos", "ref", "alt", "region", "has_sex_cols", "nonpar_xy_is_haploid", "n_xx", "n_xy"] + flag_cols + sum_cols
    return score_site_table(m[keep_cols])


def load_truth_sites(truth_path, work_dir, candidate_dfs):
    """Return set of (pos, ref, alt) truth SNPs aligned to the site tables."""
    local_truth = fetch_local(truth_path, work_dir, "truth_snp_lookup.tsv")
    ts = pd.read_csv(local_truth, sep="\t", dtype=str)
    ts["pos"] = parse_numeric_col(ts["pos"]).astype(int)

    site_keys = set()
    for df in candidate_dfs:
        site_keys.update(zip(df["pos"], df["ref"], df["alt"]))

    sample_ts = ts.drop_duplicates(["pos", "ref", "alt"]).head(20000)
    rates = {
        off: sum((p + off, r, a) in site_keys for p, r, a in zip(sample_ts["pos"], sample_ts["ref"], sample_ts["alt"]))
        for off in (0, 1)
    }
    offset = max(rates, key=rates.get)
    print(f"Truth SNP lookup offset vs site tables: {offset:+d} (matched {rates[offset]}/{len(sample_ts)} sample rows)", flush=True)
    return set(zip(ts["pos"] + offset, ts["ref"], ts["alt"]))


def build_filter_breakdown(cohorts):
    """Summarize (1) lrWGS truth x (2) ExcessHet x (3) other Step1/2 QC filters."""
    rows = []
    for cohort_label, df, _ in cohorts:
        for region in ["PAR1", "nonPAR"]:
            sub = df[df["region"] == region]
            if sub.empty:
                continue
            for (eh, in_tr, other_qc), grp in sub.groupby(["excess_het", "in_lr_truth", "f_other_qc"], dropna=False):
                rows.append({
                    "cohort": cohort_label,
                    "region": region,
                    "ExcessHet": bool(eh),
                    "in_lr_truth": bool(in_tr),
                    "fail_other_qc": bool(other_qc),
                    "n_alleles": int(len(grp)),
                    "n_fail_call_rate": int(grp["f_call_rate"].sum()),
                    "n_fail_gq": int(grp["f_gq"].sum()),
                    "n_fail_avsad": int(grp["f_avsad"].sum()),
                    "n_fail_lowqual_nohq": int(grp["f_lowqual_nohq"].sum()),
                    "n_outside_bed": int(grp["f_outside_bed"].sum()),
                })
    return pd.DataFrame(rows)


def plot_combined_3x2(df, cohort_label, n_label, region, metric, out_dir, truth_size=10.0, truth_alpha=0.60):
    """Render one 3x2 figure: rows=(XX, XY, all) x cols=(ExcessHet=False, True),
    stratifying points by (1) lrWGS truth and (3) Step1/2 other QC pass/fail."""
    sub = df[df["region"] == region].copy()
    if sub.empty:
        return None

    sex_rows = [("XX", "xx"), ("XY", "xy"), ("all", "all")]
    fig, axes = plt.subplots(3, 2, figsize=(13.5, 13.0), sharex=True, sharey=True)

    metric_titles = {
        "dev": ("Hardy-Weinberg Deviation (Observed Het − 2pq)", "HW Deviation (Observed Het − 2pq)"),
        "ratio": ("Sex-Aware Heterozygosity Ratio (H_obs / H_exp)", "Sex-Aware H_obs / H_exp"),
        "pval": ("Signed Exact HWE −log10(p) [Diploid Calls]", f"Signed −log10(p_HWE) [clipped ±{int(PVAL_CLIP)}]"),
    }
    title_prefix, ylabel = metric_titles[metric]
    fig.suptitle(f"{title_prefix} ({n_label})\nCohort: {cohort_label} | Region: {region}", fontsize=13, y=0.995)

    valid_afs = pd.concat([sub[f"af_{pfx}"] for _, pfx in sex_rows]).dropna()
    valid_afs = valid_afs[(valid_afs > 0) & (valid_afs < 1)]
    min_af = max(1e-6, float(valid_afs.min()) * 0.75) if len(valid_afs) else 1e-5
    af_range = np.logspace(np.log10(min_af), np.log10(0.999), 300)
    exp_2pq = 2.0 * af_range * (1.0 - af_range)
    min_dev = -exp_2pq
    max_ratio = 1.0 / np.maximum(af_range, 1.0 - af_range)
    xx_share = float(sub["n_xx"].iloc[0]) / float(sub["n_xx"].iloc[0] + sub["n_xy"].iloc[0])

    for r_idx, (sex_label, pfx) in enumerate(sex_rows):
        af_col = f"af_{pfx}"
        y_col = {"dev": f"hw_dev_{pfx}", "ratio": f"het_ratio_{pfx}", "pval": f"signed_p_{pfx}"}[metric]
        dip_scale = xx_share if (region == "nonPAR" and sex_label == "all") else 1.0
        max_dev = dip_scale * np.minimum(2.0 * af_range, 2.0 * (1.0 - af_range)) - exp_2pq

        for c_idx, eh_status in enumerate([False, True]):
            ax = axes[r_idx, c_idx]
            d = sub[
                (sub["excess_het"] == eh_status)
                & sub[af_col].notna()
                & (sub[af_col] > 0)
                & (sub[af_col] < 1)
                & sub[y_col].notna()
            ].copy()
            d["_y"] = d[y_col].clip(-PVAL_CLIP, PVAL_CLIP) if metric == "pval" else d[y_col]

            non_tr_fail = d[(~d["in_lr_truth"]) & d["f_other_qc"]]
            non_tr_pass = d[(~d["in_lr_truth"]) & (~d["f_other_qc"])]
            tr_fail = d[d["in_lr_truth"] & d["f_other_qc"]]
            tr_pass = d[d["in_lr_truth"] & (~d["f_other_qc"])]

            if len(non_tr_fail):
                ax.scatter(
                    non_tr_fail[af_col], non_tr_fail["_y"], alpha=0.14, s=4, rasterized=True,
                    c="lightsteelblue" if not eh_status else "lightsalmon",
                    label=f"Non-truth, fail other QC (n={len(non_tr_fail):,})",
                )
            if len(non_tr_pass):
                ax.scatter(
                    non_tr_pass[af_col], non_tr_pass["_y"], alpha=0.35, s=7, rasterized=True,
                    c="tab:blue" if not eh_status else "tab:red",
                    label=f"Non-truth, pass other QC (n={len(non_tr_pass):,})",
                )
            if len(tr_fail):
                ax.scatter(
                    tr_fail[af_col], tr_fail["_y"], alpha=truth_alpha, s=truth_size, marker="^", rasterized=True,
                    edgecolors="black", linewidths=0.25,
                    c="darkorange", label=f"Truth, fail other QC (n={len(tr_fail):,})",
                )
            if len(tr_pass):
                ax.scatter(
                    tr_pass[af_col], tr_pass["_y"], alpha=truth_alpha, s=truth_size, marker="o", rasterized=True,
                    edgecolors="black", linewidths=0.25,
                    c="gold", label=f"Truth, pass other QC (n={len(tr_pass):,})",
                )

            if metric == "dev":
                ax.axhline(0.0, color="black", linestyle="--", linewidth=1.0, label="Expected HWE (0)")
                ax.plot(af_range, max_dev, color="dimgray", linestyle=":", linewidth=1.0, label="Max diploid het")
                ax.plot(af_range, min_dev, color="dimgray", linestyle="-.", linewidth=0.9, label="Min het (0 het)")
                ax.set_ylim(-0.55, 0.55)
            elif metric == "ratio":
                ax.axhline(1.0, color="black", linestyle="--", linewidth=1.0, label="Expected HWE (1.0)")
                ax.plot(af_range, max_ratio, color="dimgray", linestyle=":", linewidth=1.0, label="Max diploid het")
                ax.axhline(0.0, color="dimgray", linestyle="-.", linewidth=0.9, label="Min het (0 het)")
                ax.set_ylim(-0.08, 2.25)
            else:
                ax.axhline(0.0, color="black", linestyle="--", linewidth=1.0, label="Expected HWE (0)")
                ax.axhline(EH_CUTOFF_NLP, color="firebrick", linestyle=":", linewidth=1.0, label="ExcessHet (p=3.4e-6)")
                ax.axhline(-EH_CUTOFF_NLP, color="steelblue", linestyle=":", linewidth=0.9, label="Het deficit (p=3.4e-6)")
                ax.set_ylim(-PVAL_CLIP * 1.08, PVAL_CLIP * 1.08)

            ax.set_xscale("log")
            ax.set_xlim(min_af, 1.02)
            ax.set_title(f"sex={sex_label} | ExcessHet={eh_status} (n={len(d):,})", fontsize=10)
            if c_idx == 0:
                ax.set_ylabel(f"[{sex_label}] {ylabel}", fontsize=9)
            if r_idx == 2:
                ax.set_xlabel("Alternate Allele Frequency (log scale)", fontsize=10)
            leg = ax.legend(loc="upper left" if metric != "dev" else "lower left", fontsize=7.0)
            for lh in getattr(leg, "legend_handles", getattr(leg, "legendHandles", [])):
                lh.set_alpha(1.0)
            ax.grid(True, which="both", alpha=0.25)

    plt.tight_layout(rect=[0, 0, 1, 0.965])
    safe_cohort = cohort_label.replace("+", "_")
    fname = f"full_cohort_hw_{metric}_{safe_cohort}_{region}.png"
    plt.savefig(os.path.join(out_dir, fname), dpi=150)
    plt.close(fig)
    return fname


def main():
    t0 = time.time()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aou-sites", required=True, help="aou.sites.tsv.bgz from export_site_tables_sex.py")
    ap.add_argument("--anvil-sites", required=True, help="anvil.sites.tsv.bgz from export_site_tables_sex.py")
    ap.add_argument("--truth-lookup", help="Optional truth_snp_lookup.tsv from vcfdist_report.py")
    ap.add_argument("--truth-size", type=float, default=10.0, help="Marker area for lrWGS truth points (default: 10)")
    ap.add_argument("--truth-alpha", type=float, default=0.60, help="Opacity for lrWGS truth points (default: 0.60)")
    ap.add_argument("--n-aou", type=int, default=414830, help="Fallback AoU N if not in TSV")
    ap.add_argument("--aou-xx-frac", type=float, default=0.606, help="Fallback AoU XX fraction if not in TSV")
    ap.add_argument("--n-anvil", type=int, default=100749, help="Fallback AnVIL N if not in TSV")
    ap.add_argument("--anvil-xx-frac", type=float, default=0.455, help="Fallback AnVIL XX fraction if not in TSV")
    ap.add_argument("--out-dir", default="hw_plots_full_cohort")
    args = ap.parse_args()

    work_dir = "full_cohort_hw_temp"
    shutil.rmtree(work_dir, ignore_errors=True)
    os.makedirs(work_dir, exist_ok=True)

    print("Loading & scoring AoU sites...", flush=True)
    aou_raw = load_site_table(args.aou_sites, work_dir, "aou.sites.tsv.bgz", args.n_aou, args.aou_xx_frac)
    aou_df = score_site_table(aou_raw)
    print(f"  AoU scored in {time.time() - t0:.1f}s", flush=True)

    print("Loading & scoring AnVIL sites...", flush=True)
    anvil_raw = load_site_table(args.anvil_sites, work_dir, "anvil.sites.tsv.bgz", args.n_anvil, args.anvil_xx_frac)
    anvil_df = score_site_table(anvil_raw)
    print(f"  AnVIL scored in {time.time() - t0:.1f}s", flush=True)

    print("Combining AoU + AnVIL sites...", flush=True)
    comb_df = combine_aou_anvil(aou_df, anvil_df)
    print(f"  Combined ({len(comb_df):,} alleles) in {time.time() - t0:.1f}s", flush=True)

    n_aou_tot = int(aou_df["n_xx"].iloc[0] + aou_df["n_xy"].iloc[0])
    n_anvil_tot = int(anvil_df["n_xx"].iloc[0] + anvil_df["n_xy"].iloc[0])
    cohorts = [
        ("AoU", aou_df, f"N={n_aou_tot:,}; XX={int(aou_df['n_xx'].iloc[0]):,}, XY={int(aou_df['n_xy'].iloc[0]):,}"),
        ("AnVIL", anvil_df, f"N={n_anvil_tot:,}; XX={int(anvil_df['n_xx'].iloc[0]):,}, XY={int(anvil_df['n_xy'].iloc[0]):,}"),
        ("AoU+AnVIL", comb_df, f"N={n_aou_tot + n_anvil_tot:,}; XX={int(comb_df['n_xx'].iloc[0]):,}, XY={int(comb_df['n_xy'].iloc[0]):,}"),
    ]

    truth_keys = (
        load_truth_sites(args.truth_lookup, work_dir, [aou_df, anvil_df, comb_df])
        if args.truth_lookup else set()
    )
    for _, df, _ in cohorts:
        df["in_lr_truth"] = [(p, r, a) in truth_keys for p, r, a in zip(df["pos"], df["ref"], df["alt"])]

    breakdown_df = build_filter_breakdown(cohorts)
    breakdown_df.to_csv(os.path.join(work_dir, "full_cohort_filter_breakdown.tsv"), sep="\t", index=False)
    print("\n=== Step1/Step2 Filter Breakdown (Truth x ExcessHet x Other QC) ===", flush=True)
    print(breakdown_df.to_string(index=False), flush=True)

    for metric in ["dev", "ratio", "pval"]:
        for cohort_label, df, n_label in cohorts:
            for region in ["PAR1", "nonPAR"]:
                fn = plot_combined_3x2(
                    df, cohort_label, n_label, region, metric, work_dir,
                    truth_size=args.truth_size, truth_alpha=args.truth_alpha,
                )
                if fn:
                    print(f"  Wrote {fn} ({time.time() - t0:.1f}s)", flush=True)

    for cohort_label, df, _ in cohorts:
        safe_cohort = cohort_label.replace("+", "_")
        df.to_csv(os.path.join(work_dir, f"full_cohort_hw_stats_{safe_cohort}.tsv.gz"), sep="\t", index=False)

    for tmp_name in ("aou.sites.tsv.bgz", "anvil.sites.tsv.bgz", "truth_snp_lookup.tsv"):
        tmp_path = os.path.join(work_dir, tmp_name)
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    if args.out_dir.startswith("gs://"):
        subprocess.run(f"gcloud storage cp {work_dir}/* {args.out_dir.rstrip('/')}/", shell=True, check=True)
    else:
        os.makedirs(args.out_dir, exist_ok=True)
        for fname in os.listdir(work_dir):
            shutil.move(os.path.join(work_dir, fname), os.path.join(args.out_dir, fname))
    shutil.rmtree(work_dir, ignore_errors=True)
    print(f"Done in {time.time() - t0:.1f}s. Outputs written to {args.out_dir}", flush=True)


if __name__ == "__main__":
    main()
