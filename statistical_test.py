#!/usr/bin/env python3
"""
Statistical test for detecting sample mislabeling and duplication
from pairwise SNV mismatch rates produced by juxtabam_run.

Tests performed:
  1. Cross-assay match test: For each biosample with paired assays,
     tests whether the observed mismatch rate is consistent with a
     true match vs. a mismatch, using leave-one-out normal distributions.
  2. Within-assay duplication test: For each pair of samples within the
     same assay type, tests whether the pair is a duplicate vs. a true
     non-duplicate.

Usage:
  python -m scripts.statistical_test \\
      --genotype_mismatch_rates output/genotype_mismatch_rates.tsv \\
      --input_info input_info.tsv

  Optional:
      --assay1 ATAC --assay2 RNA    (default: first two distinct assay types)
      --alpha 0.05                   (significance level, default 0.05)
"""

import argparse
import csv
import math
import os
import sys


# ---------------------------------------------------------------
# Math helpers (no scipy dependency)
# ---------------------------------------------------------------

def norm_pdf(x, mu, sigma):
    """Gaussian probability density function."""
    if sigma <= 0:
        return 0.0
    z = (x - mu) / sigma
    return math.exp(-0.5 * z * z) / (sigma * math.sqrt(2.0 * math.pi))


def mean_std(values):
    """Return (mean, sample_std) for a list of floats."""
    n = len(values)
    if n < 2:
        return (values[0] if n == 1 else 0.0), 0.0
    mu = sum(values) / n
    var = sum((v - mu) ** 2 for v in values) / (n - 1)
    return mu, math.sqrt(var)


# ---------------------------------------------------------------
# I/O
# ---------------------------------------------------------------

def load_input_info(path):
    """
    Read the input_info manifest (tab or comma delimited).
    Returns a list of dicts with keys: sample_id, biosample, assay, bam.
    """
    with open(path) as fh:
        sample = fh.read(4096)
        fh.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        reader = csv.DictReader(fh, dialect=dialect)
        rows = list(reader)

    required = {"sample_id", "biosample", "assay"}
    if not required.issubset(rows[0].keys()):
        raise ValueError(
            f"input_info must contain columns {required}; "
            f"found {set(rows[0].keys())}"
        )
    return rows


def load_mismatch_rates(path):
    """
    Read genotype_mismatch_rates.tsv produced by juxtabam_run.
    Returns a symmetric lookup dict: (sample1, sample2) -> percent_mismatch.
    """
    lookup = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            s1, s2 = row["sample1"], row["sample2"]
            val = float(row["percent_mismatch"])
            lookup[(s1, s2)] = val
            lookup[(s2, s1)] = val
    return lookup


def detect_assay_types(manifest_rows):
    """Return sorted list of distinct assay types from the manifest."""
    return sorted(set(r["assay"] for r in manifest_rows))


def build_sample_maps(manifest_rows):
    """
    From the input_info manifest, build:
      sid_to_bio  : sample_id -> biosample
      sid_to_assay: sample_id -> assay
      assay_bio_to_sids: (assay, biosample) -> [sample_id, ...]
    """
    sid_to_bio = {}
    sid_to_assay = {}
    assay_bio_to_sids = {}

    for r in manifest_rows:
        sid = r["sample_id"]
        bio = r["biosample"]
        assay = r["assay"]
        sid_to_bio[sid] = bio
        sid_to_assay[sid] = assay
        assay_bio_to_sids.setdefault((assay, bio), []).append(sid)

    return sid_to_bio, sid_to_assay, assay_bio_to_sids


# ---------------------------------------------------------------
# Statistical tests
# ---------------------------------------------------------------

def cross_assay_match_test(lookup, sid_to_bio, sid_to_assay, assay1, assay2):
    """
    For every pair of samples (s1, s2) where s1 is assay1 and s2 is
    assay2, classify it as a matched pair (same biosample) or mismatched
    pair (different biosample).  Each sample_id is treated individually,
    including replicates.

    For each matched pair, test whether its mismatch rate is consistent
    with the matched distribution vs. the mismatched distribution.

    Returns a list of result dicts.
    """
    # Collect sample_ids by assay
    sids_a1 = sorted(s for s, a in sid_to_assay.items() if a == assay1)
    sids_a2 = sorted(s for s, a in sid_to_assay.items() if a == assay2)

    # Partition all cross-assay pairs into matched vs mismatched
    matched_pairs = {}   # (s1, s2) -> rate
    mismatched_rates = []

    for s1 in sids_a1:
        for s2 in sids_a2:
            rate = lookup.get((s1, s2))
            if rate is None:
                continue
            if sid_to_bio[s1] == sid_to_bio[s2]:
                matched_pairs[(s1, s2)] = rate
            else:
                mismatched_rates.append(rate)

    n_matched = len(matched_pairs)
    if n_matched < 3:
        print(f"[warning] Only {n_matched} matched cross-assay pairs between "
              f"{assay1} and {assay2}; need at least 3 for meaningful statistics.")
        if n_matched < 2:
            return []

    mu_mismatch, sigma_mismatch = mean_std(mismatched_rates)
    matched_vals = list(matched_pairs.values())

    results = []
    for (s1, s2), x_i in sorted(matched_pairs.items()):
        # Leave-one-out for the matched distribution
        other_matched = [v for k, v in matched_pairs.items() if k != (s1, s2)]
        mu_match, sigma_match = mean_std(other_matched)

        f_match = norm_pdf(x_i, mu_match, sigma_match)
        f_mismatch = norm_pdf(x_i, mu_mismatch, sigma_mismatch)

        denom = f_match + f_mismatch
        p_match = f_match / denom if denom > 0 else 0.5

        bio = sid_to_bio[s1]
        label = f"{assay1}-{assay2}:{s1},{s2}" if s1 != f"{assay1}_{bio}" or s2 != f"{assay2}_{bio}" else f"{assay1}-{assay2}:{bio}"

        results.append({
            "test": "match",
            "label": label,
            "observed": x_i,
            "mu_mismatch": mu_mismatch,
            "sigma_mismatch": sigma_mismatch,
            "mu_ref": mu_match,
            "sigma_ref": sigma_match,
            "f_ref": f_match,
            "f_alt": f_mismatch,
            "p_value": p_match,
            "p_label": f"p-match:{p_match:.3g}",
        })

    return results


def within_assay_duplication_test(lookup, sid_to_bio, sid_to_assay, assay):
    """
    For every pair of samples (s1, s2) within the same assay type that
    belong to different biosamples, test whether the pair is a
    duplicate.  Each sample_id is treated individually.

    Returns a list of result dicts.
    """
    sids = sorted(s for s, a in sid_to_assay.items() if a == assay)

    # Collect all within-assay pairs from different biosamples
    pair_rates = {}
    for i, s1 in enumerate(sids):
        for s2 in sids[i + 1:]:
            if sid_to_bio[s1] == sid_to_bio[s2]:
                continue  # same biosample (replicates), skip
            rate = lookup.get((s1, s2))
            if rate is not None:
                pair_rates[(s1, s2)] = rate

    if len(pair_rates) < 3:
        print(f"[warning] Only {len(pair_rates)} cross-biosample pairs for "
              f"assay {assay}; need at least 3 for meaningful duplication "
              f"statistics.")
        if len(pair_rates) < 2:
            return []

    all_rates = list(pair_rates.values())
    mu_mismatch, sigma_mismatch = mean_std(all_rates)

    results = []
    for (s1, s2), x_ij in pair_rates.items():
        f_mismatch = norm_pdf(x_ij, mu_mismatch, sigma_mismatch)
        f_dupl = norm_pdf(x_ij, 0, sigma_mismatch)

        denom = f_mismatch + f_dupl
        p_nondupl = f_mismatch / denom if denom > 0 else 0.5

        bio1, bio2 = sid_to_bio[s1], sid_to_bio[s2]
        label = f"{assay}:{s1},{s2}" if s1 != f"{assay}_{bio1}" or s2 != f"{assay}_{bio2}" else f"{assay}:{bio1}-{bio2}"

        results.append({
            "test": "duplication",
            "label": label,
            "observed": x_ij,
            "mu_mismatch": mu_mismatch,
            "sigma_mismatch": sigma_mismatch,
            "mu_ref": 0.0,
            "sigma_ref": sigma_mismatch,
            "f_ref": f_dupl,
            "f_alt": f_mismatch,
            "p_value": p_nondupl,
            "p_label": f"p-nondupl:{p_nondupl:.3g}",
        })

    return results


# ---------------------------------------------------------------
# Output
# ---------------------------------------------------------------

def print_results(all_results, bonferroni_threshold):
    """Print sorted table of all test results."""
    all_results.sort(key=lambda r: (r["p_value"], r["observed"]))

    hdr = (f"{'p-value':>12s}  {'label':<30s}  {'observed':>8s}  "
           f"{'mu_mis':>8s}  {'sig_mis':>8s}  {'mu_ref':>8s}  "
           f"{'sig_ref':>8s}  {'f_ref':>12s}  {'f_alt':>12s}  "
           f"{'result':<30s}")
    sep = "=" * len(hdr)

    print(sep)
    print(hdr)
    print(sep)

    for r in all_results:
        flag = "***" if r["p_value"] < bonferroni_threshold else ""
        print(
            f"{r['p_value']:>12.3g}  {r['label']:<30s}  "
            f"{r['observed']:>8.3f}  "
            f"{r['mu_mismatch']:>8.3f}  {r['sigma_mismatch']:>8.3f}  "
            f"{r['mu_ref']:>8.3f}  {r['sigma_ref']:>8.3f}  "
            f"{r['f_ref']:>12.3g}  {r['f_alt']:>12.3g}  "
            f"{r['p_label']:<30s} {flag}"
        )

    flagged = [r for r in all_results if r["p_value"] < bonferroni_threshold]
    print(f"\n{sep}")
    print(f"Flagged results (p < {bonferroni_threshold:.4g}):")
    if flagged:
        for r in flagged:
            print(f"  {r['p_label']:>40s}  {r['label']}")
    else:
        print("  None")


# ---------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------

def plot_pvalue_summary(match_results, dupl_results_by_assay,
                        bonferroni_threshold, outdir):
    """
    Strip chart of all p-values on -log10 scale, grouped by test type.
    Flagged results (above -log10 threshold) are labeled.

    Parameters
    ----------
    match_results : list of result dicts from cross_assay_match_test
    dupl_results_by_assay : list of (assay_name, results_list) tuples
    bonferroni_threshold : float
    outdir : str or path
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[warning] matplotlib not available; skipping plot.")
        return

    n_groups = 1 + len(dupl_results_by_assay)
    fig, ax = plt.subplots(figsize=(2.5 + 1.8 * n_groups, 5))

    colors = ["#2ca02c", "#4393c3", "#e6550d", "#9467bd", "#8c564b"]
    threshold_neglog = -math.log10(bonferroni_threshold)

    # Group 0: match test
    x_pos = 0
    tick_labels = []
    match_label = match_results[0]["label"].split(":")[0] if match_results else "match"
    tick_labels.append(f"{match_label}\nmatch")

    match_p_neglog = [-math.log10(r["p_value"]) for r in match_results]
    ax.scatter([x_pos] * len(match_p_neglog), match_p_neglog, alpha=0.5, s=25,
               color=colors[0], zorder=3)

    for r, neglog_p in zip(match_results, match_p_neglog):
        if r["p_value"] < bonferroni_threshold:
            bio = r["label"].split(":")[-1]
            ax.annotate(bio, (x_pos, neglog_p),
                        fontsize=9, ha="left", color="#d62728",
                        xytext=(8, 0), textcoords="offset points")

    # Groups 1..N: duplication tests per assay
    for idx, (assay_name, dupl_results) in enumerate(dupl_results_by_assay):
        x_pos = idx + 1
        tick_labels.append(f"{assay_name}\nduplication")
        c = colors[(idx + 1) % len(colors)]

        dupl_p_neglog = [-math.log10(r["p_value"]) for r in dupl_results]
        ax.scatter([x_pos] * len(dupl_p_neglog), dupl_p_neglog, alpha=0.35, s=18,
                   color=c, zorder=3)

        for r, neglog_p in zip(dupl_results, dupl_p_neglog):
            if r["p_value"] < bonferroni_threshold:
                pair = r["label"].split(":")[-1]
                ax.annotate(pair, (x_pos, neglog_p),
                            fontsize=9, ha="left", color="#d62728",
                            xytext=(8, 0), textcoords="offset points")

    ax.axhline(threshold_neglog, color="black", linestyle="--",
               linewidth=1, label=f"Bonferroni",
               zorder=2)

    ax.set_xticks(range(len(tick_labels)))
    ax.set_xticklabels(tick_labels, fontsize=10)
    ax.set_ylabel("-log10(p-value)")
    ax.set_title("Statistical test results")
    ax.legend(fontsize=9, loc="upper right")
    ax.set_xlim(-0.5, len(tick_labels) - 0.5)
    plt.tight_layout()

    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, "pvalue_summary.png")
    plt.savefig(out_path, dpi=200)
    plt.close()
    print(f"\nSaved plot: {out_path}")


# ---------------------------------------------------------------
# CLI
# ---------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Detect sample mislabeling and duplication from pairwise "
            "SNV mismatch rates produced by juxtabam_run."
        )
    )
    parser.add_argument(
        "--genotype_mismatch_rates", required=True,
        help="Path to genotype_mismatch_rates.tsv from juxtabam_run."
    )
    parser.add_argument(
        "--input_info", required=True,
        help="Path to input_info manifest (same file used for juxtabam_run)."
    )
    parser.add_argument(
        "--assay1", default=None,
        help="First assay type for cross-assay match test "
             "(default: first assay type found in input_info)."
    )
    parser.add_argument(
        "--assay2", default=None,
        help="Second assay type for cross-assay match test "
             "(default: second assay type found in input_info)."
    )
    parser.add_argument(
        "--alpha", type=float, default=0.05,
        help="Significance level for Bonferroni correction (default: 0.05)."
    )
    parser.add_argument(
        "--outdir", default=None,
        help="Output directory for plots. If not specified, no plot is generated."
    )

    return parser.parse_args()


def main():
    args = parse_args()

    # Load data
    manifest_rows = load_input_info(args.input_info)
    lookup = load_mismatch_rates(args.genotype_mismatch_rates)
    sid_to_bio, sid_to_assay, assay_bio_to_sids = build_sample_maps(manifest_rows)

    # Resolve assay types
    all_assays = detect_assay_types(manifest_rows)
    print(f"Assay types found in input_info: {all_assays}")

    if len(all_assays) < 2:
        print("[error] Need at least two distinct assay types for "
              "cross-assay match testing.")
        sys.exit(1)

    assay1 = args.assay1 if args.assay1 else all_assays[0]
    assay2 = args.assay2 if args.assay2 else all_assays[1]

    if assay1 not in all_assays:
        print(f"[error] --assay1 '{assay1}' not found in input_info. "
              f"Available: {all_assays}")
        sys.exit(1)
    if assay2 not in all_assays:
        print(f"[error] --assay2 '{assay2}' not found in input_info. "
              f"Available: {all_assays}")
        sys.exit(1)
    if assay1 == assay2:
        print("[error] --assay1 and --assay2 must be different assay types.")
        sys.exit(1)

    # Count shared biosamples
    bios_a1 = {bio for (a, bio) in assay_bio_to_sids if a == assay1}
    bios_a2 = {bio for (a, bio) in assay_bio_to_sids if a == assay2}
    shared = sorted(bios_a1 & bios_a2)
    print(f"Testing {assay1} vs {assay2}: "
          f"{len(shared)} shared biosamples\n")

    # Run tests
    match_results = cross_assay_match_test(
        lookup, sid_to_bio, sid_to_assay, assay1, assay2
    )
    dupl_results_a1 = within_assay_duplication_test(
        lookup, sid_to_bio, sid_to_assay, assay1
    )
    dupl_results_a2 = within_assay_duplication_test(
        lookup, sid_to_bio, sid_to_assay, assay2
    )

    all_results = match_results + dupl_results_a1 + dupl_results_a2

    n_tests = len(all_results)
    bonferroni_threshold = args.alpha / n_tests if n_tests > 0 else args.alpha
    print(f"Total tests: {n_tests} "
          f"({len(match_results)} match + "
          f"{len(dupl_results_a1)} {assay1}-dupl + "
          f"{len(dupl_results_a2)} {assay2}-dupl)")
    print(f"Bonferroni threshold ({args.alpha}/{n_tests}): "
          f"{bonferroni_threshold:.4g}\n")

    print_results(all_results, bonferroni_threshold)

    if args.outdir:
        plot_pvalue_summary(
            match_results,
            [(assay1, dupl_results_a1), (assay2, dupl_results_a2)],
            bonferroni_threshold,
            args.outdir,
        )


if __name__ == "__main__":
    main()
