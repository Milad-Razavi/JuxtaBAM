"""
Statistical tests for sample mislabeling and duplication, from the pairwise SNV
genotype comparisons produced by juxtabam_run.

For each pair of samples the pipeline reports n SNV sites and m discordant (mismatched) genotypes.
Discordance at a site is a Bernoulli event, but the underlying
discordance probability varies from pair to pair and assay to assay.
Letting the per-pair probability be random gives the beta-binomial:

    m ~ BetaBinomial(n, alpha, beta)

alpha and beta are fit by maximum likelihood to the other pairs of the same
family, excluding the pair under test (leave one out).

Two tests, each one-sided because the alternative can move the mismatch rate in
only one direction:

  Cross-assay mismatch test, for every unordered pair of assay types.
    Null: the labeled pair is the same biosample. A swap raises the mismatch.
    p = P(M >= m), the upper tail.

  Within-assay duplication test, for every assay type.
    Null: the pair is two different biosamples. A duplicate lowers the mismatch.
    p = P(M <= m), the lower tail.

Bonferroni correction is applied within each family of tests.

Usage:
  python statistical_test.py \\
      --genotype_mismatch_rates output/genotype_mismatch_rates.tsv \\
      --input_info input_info.tsv \\
      --outdir stats_tables/ [--assays ATAC,RNA] [--out results.tsv]
"""

import argparse
import os
import sys
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import betabinom

MIN_REF_PAIRS = 3
LOG_BOUNDS = ((-12.0, 25.0), (-12.0, 25.0))
P_UNDERFLOW = 1e-300


# ---------------------------------------------------------------
# Beta-binomial null
# ---------------------------------------------------------------

def neg_log_lik(log_params, m, n):
    alpha, beta = np.exp(log_params)
    return -betabinom.logpmf(m, n, alpha, beta).sum()


def fit_null(m, n):
    """
    Maximum likelihood fit of (alpha, beta) to reference counts.

    Returns (None, None) only if the optimizer returns a non-finite point.
    """
    res = minimize(neg_log_lik, (0.0, 0.0), args=(m, n),
                   bounds=LOG_BOUNDS, method="L-BFGS-B")
    if not np.all(np.isfinite(res.x)) or not np.isfinite(res.fun):
        return None, None
    alpha, beta = np.exp(res.x)
    return float(alpha), float(beta)


def tail_pvalue(m, n, alpha, beta, upper):
    """
    P(M >= m | n) if upper, else P(M <= m | n).

    scipy's sf() is P(M > k).
    """
    if not upper:
        return float(betabinom.cdf(m, n, alpha, beta))
    if m < 1:
        return 1.0
    return float(betabinom.sf(m - 1, n, alpha, beta))


def format_p(p):
    """Render a p-value in scientific notation, or < 1e-300 if underflow."""
    return f"{p:.4g}" if p > P_UNDERFLOW else f"< {P_UNDERFLOW:.0e}"


# ---------------------------------------------------------------
# I/O
# ---------------------------------------------------------------

def load_inputs(rates_path, info_path):
    """Read the manifest and the pairwise counts, annotated with biosample and assay."""
    info = pd.read_csv(info_path, sep=None, engine="python")
    missing = {"sample_id", "biosample", "assay"} - set(info.columns)
    if missing:
        raise ValueError(f"input_info is missing columns: {sorted(missing)}")

    pairs = pd.read_csv(rates_path, sep=None, engine="python")
    pairs.columns = [c.strip() for c in pairs.columns]
    need = {"sample1", "sample2", "total_sites", "mismatches"}
    if not need.issubset(pairs.columns):
        raise ValueError(
            f"genotype_mismatch_rates must contain {sorted(need)}; found "
            f"{sorted(pairs.columns)}."
        )
    if "percent_mismatch" not in pairs.columns:
        pairs["percent_mismatch"] = 100 * pairs.mismatches / pairs.total_sites

    meta = info.set_index("sample_id")[["biosample", "assay"]]
    for i in ("1", "2"):
        pairs = pairs.join(meta.add_suffix(i), on=f"sample{i}")
    return pairs.dropna(subset=["biosample1", "biosample2"])


# ---------------------------------------------------------------
# Tests
# ---------------------------------------------------------------

def run_test(df, upper, family, test_type, min_sites):
    """Leave-one-out beta-binomial tail test over the rows of df."""
    if len(df) <= MIN_REF_PAIRS:
        print(f"[warning] Family '{family}' has {len(df)} pair(s); more than "
              f"{MIN_REF_PAIRS} are needed to fit a null. Skipping this family.")
        return []

    m = df.mismatches.to_numpy(dtype=int)
    n = df.total_sites.to_numpy(dtype=int)

    out = []
    for i, row in enumerate(df.itertuples()):
        alpha, beta = fit_null(np.delete(m, i), np.delete(n, i))
        if alpha is None:
            print(f"[warning] Null fit returned a non-finite point for '{family}' "
                  f"excluding {row.sample1},{row.sample2}. Skipping this pair.")
            continue
        out.append(dict(
            family=family,
            test_type=test_type,
            sample1=row.sample1,
            sample2=row.sample2,
            biosample1=row.biosample1,
            biosample2=row.biosample2,
            obs_mismatch_percent=row.percent_mismatch,
            mismatches=m[i],
            total_sites=n[i],
            alpha=alpha,
            beta=beta,
            exp_mismatch_percent=100 * alpha / (alpha + beta),
            p_value=tail_pvalue(m[i], n[i], alpha, beta, upper),
            low_sites=n[i] < min_sites,
        ))
    return out


def cross_assay_mismatch_test(pairs, a1, a2, min_sites):
    """Same biosample, one sample of each assay type. Upper tail."""
    sel = pairs[
        (pairs.biosample1 == pairs.biosample2)
        & (((pairs.assay1 == a1) & (pairs.assay2 == a2))
           | ((pairs.assay1 == a2) & (pairs.assay2 == a1)))
    ]
    return run_test(sel, True, f"{a1}__{a2}", "cross_assay_mismatch", min_sites)


def within_assay_duplication_test(pairs, assay, min_sites):
    """Same assay type, different biosamples. Lower tail."""
    sel = pairs[
        (pairs.assay1 == assay) & (pairs.assay2 == assay)
        & (pairs.biosample1 != pairs.biosample2)
    ]
    return run_test(sel, False, assay, "within_assay_duplication", min_sites)


def run_all_tests(pairs, assays, min_sites):
    results = []
    for a1, a2 in combinations(assays, 2):
        results += cross_assay_mismatch_test(pairs, a1, a2, min_sites)
    for assay in assays:
        results += within_assay_duplication_test(pairs, assay, min_sites)
    return results


def apply_bonferroni(results, alpha):
    """Bonferroni within each family: threshold = alpha / (tests in that family)."""
    df = pd.DataFrame(results)
    df["n_tests"] = df.groupby(["test_type", "family"])["family"].transform("size")
    df["threshold"] = alpha / df.n_tests
    df["significant"] = df.p_value < df.threshold
    df["adjusted_p"] = np.minimum(1.0, df.p_value * df.n_tests)
    return df.sort_values(["test_type", "family", "p_value"]).reset_index(drop=True)


# ---------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------

TEST_LEGEND = {
    "cross_assay_mismatch": [
        "Cross-assay mismatch test: null is 'same biosample'; upper tail P(M >= m).",
        "A small p means the pair is more discordant than a true match should be,",
        "suggesting a swap or mislabel.",
    ],
    "within_assay_duplication": [
        "Within-assay duplication test: null is 'different biosamples'; lower tail",
        "P(M <= m). A small p means the pair is more concordant than two distinct",
        "sources should be, suggesting a duplicate.",
    ],
}


def _header(alpha, assays):
    return [
        "=" * 100,
        "JuxtaBAM statistical tests: sample mislabeling and duplication",
        "=" * 100,
        "",
        f"Assay types        : {', '.join(assays)}",
        f"Significance level : alpha = {alpha}",
        "Correction         : Bonferroni, applied within each test family",
        "",
    ]


def family_table(test_type, family, sub, alpha):
    n_tests = int(sub.n_tests.iloc[0])
    thr = float(sub.threshold.iloc[0])
    lines = [
        f"Family: {test_type} / {family}",
        f"  tests in family      : {n_tests}",
        f"  Bonferroni threshold : {alpha} / {n_tests} = {thr:.6g}",
        "",
    ]
    hdr = (f"  {'p_value':>12s}  {'adjusted_p':>12s}  {'sample1':<32s} {'sample2':<32s} "
           f"{'obs_mismatch_percent':>20s} {'mism':>9s} {'total_sites':>11s} "
           f"{'exp_mismatch_percent':>20s} {'alpha':>10s} {'beta':>12s}  flag")
    lines.append(hdr)
    lines.append("  " + "-" * (len(hdr) - 2))
    for r in sub.itertuples():
        flag = "***" if r.significant else ("low_n" if r.low_sites else "")
        lines.append(
            f"  {format_p(r.p_value):>12s}  {format_p(r.adjusted_p):>12s}  "
            f"{r.sample1:<32s} {r.sample2:<32s} "
            f"{r.obs_mismatch_percent:>20.3f} {r.mismatches:>9d} {r.total_sites:>11d} "
            f"{r.exp_mismatch_percent:>20.3f} {r.alpha:>10.4f} {r.beta:>12.4f}  {flag}"
        )
    lines.append("")
    lines.append(f"  *** marks p below the family Bonferroni threshold ({thr:.6g}).")
    if sub.low_sites.any():
        lines.append("  low_n marks a comparison resting on few usable SNV sites; its null is")
        lines.append("  wide, so a non-significant result there is weak evidence.")
    lines.append("")
    return lines


def summary_table(df, alpha, assays):
    lines = _header(alpha, assays)
    lines += ["=" * 100, "FAMILIES TESTED", "=" * 100, ""]

    # Size the text columns to the actual data, so long assay or family names
    # never overflow a fixed width and shift the columns after them.
    groups = list(df.groupby(["test_type", "family"], sort=False))
    tt_w = max([len("test type")] + [len(tt) for (tt, _), _ in groups])
    fam_w = max([len("family")] + [len(fam) for (_, fam), _ in groups])

    lines.append(f"  {'test type':<{tt_w}s}  {'family':<{fam_w}s}  {'tests':>7s}  "
                 f"{'threshold':>12s}  {'flagged':>8s}")
    lines.append("  " + "-" * (tt_w + fam_w + 7 + 12 + 8 + 8))
    for (test_type, family), sub in groups:
        lines.append(f"  {test_type:<{tt_w}s}  {family:<{fam_w}s}  "
                     f"{int(sub.n_tests.iloc[0]):>7d}  "
                     f"{float(sub.threshold.iloc[0]):>12.6g}  "
                     f"{int(sub.significant.sum()):>8d}")
    lines += ["", "=" * 100, "FLAGGED COMPARISONS", "=" * 100, ""]

    hits = df[df.significant]
    if hits.empty:
        lines.append("  None. No comparison falls below its family Bonferroni threshold.")
        lines.append("")
    for r in hits.itertuples():
        reason = ("elevated mismatch, possible swap or mislabel"
                  if r.test_type == "cross_assay_mismatch"
                  else "unexpectedly low mismatch, possible duplicate")
        lines.append(f"  [{r.test_type} / {r.family}]  {r.sample1} vs {r.sample2}")
        lines.append(f"      obs_mismatch_percent {r.obs_mismatch_percent:.3f}%  "
                     f"(exp_mismatch_percent {r.exp_mismatch_percent:.3f}%),  "
                     f"{r.mismatches} / {r.total_sites} total_sites,  p = {format_p(r.p_value)}")
        lines.append(f"      {reason}")
        lines.append("")

    weak = df[(~df.significant) & df.low_sites]
    if not weak.empty:
        lines.append(f"[note] {len(weak)} comparison(s) rest on fewer than the minimum "
                     f"site count; their")
        lines.append("       nulls are wide, so a non-significant result there is weak "
                     "evidence.")
        lines.append("")
    return lines


def write_tables(df, alpha, assays, outdir):
    os.makedirs(outdir, exist_ok=True)
    written = []

    path = os.path.join(outdir, "summary.txt")
    with open(path, "w") as fh:
        fh.write("\n".join(summary_table(df, alpha, assays)) + "\n")
    written.append(path)

    for (test_type, family), sub in df.groupby(["test_type", "family"], sort=False):
        safe = family.replace("/", "_")
        path = os.path.join(outdir, f"{test_type}__{safe}.txt")
        lines = (_header(alpha, assays) + TEST_LEGEND[test_type] + [""]
                 + family_table(test_type, family, sub, alpha))
        with open(path, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        written.append(path)
    return written


def main():
    ap = argparse.ArgumentParser(
        description="Detect sample mislabeling and duplication from JuxtaBAM "
                    "pairwise SNV genotype comparisons.")
    ap.add_argument("--genotype_mismatch_rates", required=True)
    ap.add_argument("--input_info", required=True)
    ap.add_argument("--outdir", required=True,
                    help="Directory for the per-family tables and the summary.")
    ap.add_argument("--assays", default=None,
                    help="Comma-separated assay types to test. Default: all found.")
    ap.add_argument("--alpha", type=float, default=0.05,
                    help="Family-wise significance level (default 0.05).")
    ap.add_argument("--min_sites", type=int, default=50,
                    help="Annotate comparisons resting on fewer usable SNV sites.")
    ap.add_argument("--out", default=None,
                    help="Also write the combined result table as TSV.")
    args = ap.parse_args()

    pairs = load_inputs(args.genotype_mismatch_rates, args.input_info)

    found = sorted(pd.unique(pd.concat([pairs.assay1, pairs.assay2]).dropna()))
    if args.assays:
        assays = [a.strip() for a in args.assays.split(",") if a.strip()]
        unknown = set(assays) - set(found)
        if unknown:
            sys.exit(f"[error] assay type(s) not in the data: {sorted(unknown)}. "
                     f"Available: {found}")
    else:
        assays = found

    print(f"Assay types found : {found}")
    print(f"Assay types tested: {assays}")
    if len(assays) < 2:
        print("[warning] Fewer than two assay types; no cross-assay test is possible.")

    results = run_all_tests(pairs, assays, args.min_sites)
    if not results:
        sys.exit("[error] No testable pairs found.")

    df = apply_bonferroni(results, args.alpha)

    print()
    print("\n".join(summary_table(df, args.alpha, assays)))

    for path in write_tables(df, args.alpha, assays, args.outdir):
        print(f"Wrote: {path}")
    if args.out:
        df.to_csv(args.out, sep="\t", index=False)
        print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()