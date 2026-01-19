import os
import subprocess
import pandas as pd
from pathlib import Path
from itertools import combinations
import time
from juxtabam.utils import run, build_regions_arg


def slurm_wait(job_ids, poll_seconds=30):
    """
    Wait until all job IDs are no longer present in squeue.
    """
    if not job_ids:
        return

    job_ids = [str(j).strip() for j in job_ids if str(j).strip()]
    job_ids = [j.split(";")[0] for j in job_ids]  # sbatch --parsable safety
    job_ids_str = ",".join(job_ids)

    print(f"[info] Waiting on Slurm jobs: {job_ids_str}")

    while True:
        # squeue exits 0 even if none of the jobs exist; output will be empty
        cmd = f"squeue -h -j {job_ids_str} -o %i"
        out = subprocess.check_output(cmd, shell=True).decode().strip().splitlines()
        active = set([x.strip() for x in out if x.strip()])

        if not active:
            print("[info] All Slurm jobs finished.")
            return

        print(f"[info] Still running: {len(active)} jobs. Sleeping {poll_seconds}s.")
        time.sleep(poll_seconds)


def slurm_submit(script_path):
    """
    Submit an sbatch script and return the job ID.
    Requires sbatch to be available in PATH.
    """
    cmd = f"sbatch --parsable {script_path}"
    out = subprocess.check_output(cmd, shell=True).decode().strip()
    # --parsable can return "jobid" or "jobid;cluster"
    job_id = out.split(";")[0]
    return job_id



def hc_slurm_worker(args):
    """
    Slurm worker mode: run HaplotypeCaller for exactly one sample_id, then exit.
    """
    gatk = args.gatk_path

    outdir = Path(args.output_dir)
    outdir.mkdir(exist_ok=True, parents=True)

    if args.tmp_dir is None:
        tmp_dir = outdir / "tmp"
    else:
        tmp_dir = Path(args.tmp_dir)
    tmp_dir.mkdir(exist_ok=True, parents=True)

    gatk_java_opts = f'{args.java_mem} -Djava.io.tmpdir={tmp_dir}'

    # Regions
    regions_arg = build_regions_arg(args.regions)
    # if args.regions:
    #     regions = parse_regions(args.regions)
    #     if len(regions) == 1 and os.path.exists(regions[0]):
    #         regions_arg = f"-L {regions[0]}"
    #     else:
    #         regions_arg = " ".join([f"-L {r}" for r in regions])

    gvcf_dir = outdir / "gvcf"
    gvcf_dir.mkdir(exist_ok=True, parents=True)

    # Load processed manifest
    manifest_path = Path(args.processed_manifest)
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"Processed BAM manifest not found: {manifest_path}\n"
            f"This file is created by the main run before Slurm submission."
        )

    dfp = pd.read_csv(manifest_path, sep="\t")
    if "sample_id" not in dfp.columns or "processed_bam" not in dfp.columns:
        raise ValueError(f"Processed manifest must contain: sample_id, processed_bam. Found: {list(dfp.columns)}")

    sub = dfp[dfp["sample_id"] == args.slurm_sample_id]
    if sub.empty:
        raise ValueError(f"Sample id not found in processed manifest: {args.slurm_sample_id}")

    bam = str(sub["processed_bam"].iloc[0])
    sid = args.slurm_sample_id
    gvcf = gvcf_dir / f"{sid}.g.vcf.gz"

    if gvcf.exists():
        print(f"[warning] GVCF already exists for {sid}: {gvcf}")
        print("          Skipping HaplotypeCaller for this sample.")
        return

    cmd = f"""
{gatk} --java-options "{gatk_java_opts}" HaplotypeCaller \
  -R {args.reference} \
  -I {bam} \
  --emit-ref-confidence GVCF \
  -O {gvcf} {regions_arg} \
  --native-pair-hmm-threads {args.threads}
"""
    run(cmd)


def pairwise_bcftools_slurm_worker(args):
    """
    Slurm worker: run bcftools query for exactly one sample pair.
    Writes raw output to a tmp file.
    """
    bcftools = args.bcftools_path
    out_path = Path(args._pair_out)

    s1, s2 = args._pair.split(",")

    out_path.parent.mkdir(exist_ok=True, parents=True)

    cmd = f"""
{bcftools} query \
  -s {s1},{s2} \
  -f '%CHROM\\t%POS\\t[%GT\\t%DP\\t]\\n' \
  {args.pass_vcf} > {out_path}
"""
    run(cmd)


def write_processed_manifest(df, path):
    """
    Save processed_bam mapping so Slurm worker can run HaplotypeCaller without redoing steps.
    """
    cols = ["sample_id", "processed_bam"]
    df[cols].to_csv(path, sep="\t", index=False)