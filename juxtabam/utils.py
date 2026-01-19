import subprocess
import os

def run(cmd):
    print("\n>>>", cmd, "\n")
    subprocess.check_call(cmd, shell=True)
    

def bam_is_valid(bam, samtools):
    """
    Return True if BAM is readable & complete.
    """
    try:
        subprocess.check_call(
            f"{samtools} quickcheck -q {bam}",
            shell=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        return True
    except subprocess.CalledProcessError:
        return False


def index_is_valid(bam):
    """
    Check if the BAM index exists.
    """
    return os.path.exists(bam + ".bai")


def bam_has_good_rg(bam, samtools):
    """
    Return True if BAM contains at least one @RG line.
    """
    try:
        out = subprocess.check_output(
            f"{samtools} view -H {bam}",
            shell=True
        ).decode()

        for line in out.splitlines():
            if line.startswith("@RG"):
                return True
        return False

    except Exception:
        return False


def parse_regions(region_string):
    if region_string is None:
        return None
    regions = [r.strip() for r in region_string.split(",")]
    return regions


def build_regions_arg(regions, exists_fn=os.path.exists):
    regions_arg = ""
    if regions:
        parsed = parse_regions(regions)
        if len(parsed) == 1 and exists_fn(parsed[0]):
            regions_arg = f"-L {parsed[0]}"
        else:
            regions_arg = " ".join([f"-L {r}" for r in parsed])
    return regions_arg


def normalize_gt(gt):
    """
    Normalize genotype representation so that:
    0/1 == 1/0 == 0|1 == 1|0
    """
    if gt in (".", "./.", ".|."):
        return None

    gt = gt.replace("|", "/")
    alleles = gt.split("/")

    if len(alleles) != 2:
        return None
    if "." in alleles:
        return None

    alleles = sorted(alleles)
    return "/".join(alleles)

