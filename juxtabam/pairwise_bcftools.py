# pairwise_bcftools.py
import subprocess
from juxtabam.utils import normalize_gt


def build_bcftools_pairwise_query(bcftools, pass_vcf, s1, s2):
    return (
        f"{bcftools} query "
        f"-s {s1},{s2} "
        f"-f '%CHROM\\t%POS\\t[%GT\\t%DP\\t]\\n' "
        f"{pass_vcf}"
    )


def run_bcftools_pairwise_query(cmd):
    return subprocess.check_output(cmd, shell=True).decode().splitlines()


def parse_bcftools_pairwise_output(
    lines,
    min_dp,
    min_sites_warning,
    s1,
    s2,
    warn_fn=print,
):
    total = 0
    mismatches = 0

    for line in lines:
        fields = line.split("\t")
        if len(fields) < 6:
            continue

        gt1, dp1 = fields[2], fields[3]
        gt2, dp2 = fields[4], fields[5]

        ngt1 = normalize_gt(gt1)
        ngt2 = normalize_gt(gt2)
        if ngt1 is None or ngt2 is None:
            continue

        try:
            dp1 = int(dp1)
            dp2 = int(dp2)
        except ValueError:
            continue

        if dp1 < min_dp or dp2 < min_dp:
            continue

        total += 1
        if ngt1 != ngt2:
            mismatches += 1

    if total < min_sites_warning:
        warn_fn(
            f"[warning] Only {total} usable SNPs for comparison {s1} vs {s2}. Interpret cautiously."
        )

    percent = (100 * mismatches / total) if total > 0 else None
    return total, mismatches, percent
