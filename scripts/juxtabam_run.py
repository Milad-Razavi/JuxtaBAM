#!/usr/bin/env python3
import argparse
import os
import subprocess
from pathlib import Path
from itertools import combinations
import datetime
import platform
import time

import pandas as pd

from juxtabam.pairwise_bcftools import (
    build_bcftools_pairwise_query,
    run_bcftools_pairwise_query,
    parse_bcftools_pairwise_output,
)
from juxtabam.slurm_utils import (
    slurm_wait,
    slurm_submit,
    hc_slurm_worker,
    pairwise_bcftools_slurm_worker,
    write_processed_manifest,
)
from juxtabam.utils import (
    run,
    index_is_valid,
    bam_is_valid,
    bam_has_good_rg,
    build_regions_arg
)


MIN_SITES_WARNING = 100
RUN_INFO_FILE = "run_info.txt"
PROCESSED_MANIFEST_FILE = "processed_bams.tsv"
PAIRWISE_TMP_DIR = "pairwise_tmp"


def load_and_validate_input(args, samtools):
    df = pd.read_csv(
        args.input_info,
        sep=None,
        engine="python"
    )

    assert set(["sample_id", "biosample", "assay", "bam"]).issubset(df.columns)
    
    if df.shape[0] < 2:
        raise ValueError(
            f"Input file must contain at least two samples. Found only {df.shape[0]}."
        )

    ########################
    # Basic sanity checks
    ########################
    if df["sample_id"].duplicated().any():
        dups = df[df["sample_id"].duplicated()]["sample_id"].tolist()
        raise ValueError(f"Duplicate sample_id entries found: {dups}")

    # BAM validity checks
    for bam in df["bam"]:
        if not os.path.exists(bam):
            raise FileNotFoundError(bam)

        if not index_is_valid(bam):
            raise FileNotFoundError(f"Missing BAM index: {bam}.bai")

        if not bam_is_valid(bam, samtools):
            raise RuntimeError(
                f"BAM appears corrupted or truncated: {bam}\n"
                f"(samtools quickcheck failed)"
            )

        if not bam_has_good_rg(bam, samtools):
            print(f"[warning] No @RG line found in BAM header: {bam}")
            print("          AddOrReplaceReadGroups will fix this.")

    # Reference integrity
    fai = args.reference + ".fai"
    dictfile = args.reference.replace(".fa", "").replace(".fasta", "") + ".dict"
    if not os.path.exists(fai):
        raise FileNotFoundError(f"Missing FASTA index: {fai}")
    if not os.path.exists(dictfile):
        raise FileNotFoundError(f"Missing sequence dictionary: {dictfile}")

    ########################
    # Assay imbalance info
    ########################
    biosample_counts = df.groupby("biosample")["assay"].nunique()
    missing = biosample_counts[biosample_counts < 2]

    if len(missing):
        print("[info] Some biosamples may not have multiple assays available:")
        for b in missing.index:
            print(f"  {b}")
    return df


def prepare_output_dirs(args):
    outdir = Path(args.output_dir)
    outdir.mkdir(exist_ok=True, parents=True)

    if args.tmp_dir is None:
        tmp_dir = outdir / "tmp"
    else:
        tmp_dir = Path(args.tmp_dir)

    tmp_dir.mkdir(exist_ok=True, parents=True)

    print(f"[info] Using GATK tmp dir: {tmp_dir}")

    gvcf_dir = outdir / "gvcf"
    db_dir = outdir / "gendb"
    joint_dir = outdir / "joint"
    rg_bam_dir = outdir / "bam_rg"
    split_bam_dir = outdir / "bam_split"

    for d in [gvcf_dir, joint_dir, rg_bam_dir, split_bam_dir]:
        d.mkdir(exist_ok=True, parents=True)

    existing = [p for p in gvcf_dir.glob("*.g.vcf.gz")]
    if existing:
        print(f"[warning] {len(existing)} GVCF files already exist in {gvcf_dir}")
        print("          They will be reused unless deleted.")

    return outdir, tmp_dir, gvcf_dir, db_dir, joint_dir, rg_bam_dir, split_bam_dir


def write_run_info_start(outdir, args, split_assays, start_timestamp):
    with open(outdir / RUN_INFO_FILE, "w") as f:
        f.write(f"Start time: {start_timestamp}\n")
        f.write(f"Command: {' '.join(os.sys.argv)}\n")
        f.write(f"Reference: {args.reference}\n")
        f.write(f"Regions: {args.regions}\n")
        f.write(f"min_dp: {args.min_dp}\n")
        f.write(f"QD >= {args.min_qd}\n")
        f.write(f"QUAL >= {args.min_qual}\n")
        f.write(f"Threads: {args.threads}\n")
        f.write(f"Java memory: {args.java_mem}\n")
        f.write(f"Platform (RGPL): {args.platform}\n")
        f.write(f"SplitNCigar assays: {','.join(split_assays)}\n")
        f.write(f"Python: {platform.python_version()}\n")


def write_run_info_end(outdir, start_time):
    end_time = time.time()
    end_timestamp = datetime.datetime.now()
    runtime_seconds = int(end_time - start_time)
    
    hours = runtime_seconds // 3600
    minutes = (runtime_seconds % 3600) // 60
    seconds = runtime_seconds % 60

    with open(outdir / RUN_INFO_FILE, "a") as f:
        f.write(f"Finish time: {end_timestamp}\n")
        f.write(f"Runtime (seconds): {runtime_seconds}\n")
        f.write(f"Runtime (h:m:s): {hours}:{minutes:02d}:{seconds:02d}\n")


def process_bams(
    df,
    args,
    tools,
    gatk_java_opts,
    split_assays,
    rg_bam_dir,
    split_bam_dir,
    joint_dir,
    regions_arg,
):
    df["processed_bam"] = None

    for idx, row in df.iterrows():
        sid = row["sample_id"]
        assay = str(row["assay"])
        raw_bam = row["bam"]

        rg_bam = rg_bam_dir / f"{sid}.rg.bam"
        split_bam = split_bam_dir / f"{sid}.split.bam"

        # Step 0a: AddOrReplaceReadGroups
        if args.skip_AddOrReplaceRG:
            print(f"[info] Skipping AddOrReplaceReadGroups for {sid}")
            rg_bam = Path(raw_bam)  # use original BAM
        else:
            if rg_bam.exists():
                print(f"[warning] RG BAM already exists for {sid}: {rg_bam}")
                print("          Skipping AddOrReplaceReadGroups for this sample.")
            else:
                cmd = f"""
{tools["gatk"]} --java-options "{gatk_java_opts}" AddOrReplaceReadGroups \
  -I {raw_bam} \
  -O {rg_bam} \
  --RGID rgid_{sid} \
  --RGLB lib_{sid} \
  --RGPL {args.platform} \
  --RGPU unit_{sid} \
  --RGSM {sid} \
  --VALIDATION_STRINGENCY LENIENT
"""
                run(cmd)
                run(f"{tools['samtools']} index {rg_bam}")

        # Step 0b: SplitNCigarReads for selected assays
        assay_lower = assay.lower()

        if args.skip_SplitNCigar:
            print(f"[info] Skipping SplitNCigarReads for {sid}")
            final_bam = rg_bam
        else:
            if assay_lower in split_assays:
                if split_bam.exists():
                    print(f"[warning] SplitNCigar BAM already exists for {sid}: {split_bam}")
                    print("          Reusing existing SplitNCigarReads output.")
                else:
                    cmd = f"""
{tools["gatk"]} --java-options "{gatk_java_opts}" SplitNCigarReads \
  -R {args.reference} \
  -I {rg_bam} \
  -O {split_bam} {regions_arg}
"""
                    run(cmd)
                    print(f"[info] Running SplitNCigarReads on sample {sid} with regions: {args.regions}")
                    run(f"{tools['samtools']} index {split_bam}")

                final_bam = split_bam
            else:
                final_bam = rg_bam

        df.loc[df.sample_id == sid, "processed_bam"] = str(final_bam)

    # Save processed bam mapping for Slurm HaplotypeCaller worker jobs
    processed_manifest = joint_dir / PROCESSED_MANIFEST_FILE
    write_processed_manifest(df, processed_manifest)
    return df, processed_manifest


def parse_args():
    
    parser = argparse.ArgumentParser(
        description="Cross-assay identity checking using GATK joint SNP genotyping"
    )

    # -------------------------
    # Main arguments
    # -------------------------
    main = parser.add_argument_group("Main arguments")

    main.add_argument("--input_info", required=True,
        help="Tab-separated file containing header line: sample_id biosample assay bam")

    main.add_argument("--reference", required=True,
        help="Reference genome FASTA (.fa) file. Must be indexed by [samtools faidx] and [gatk CreateSequenceDictionary]. The 3 files must be in the same directory. Example: --reference hg38.fa while hg38.fa.fai and hg38.dict are in the same directory as hg38.fa")

    main.add_argument("--output_dir", required=True,
        help="Base output directory")

    main.add_argument("--regions", default="chr19",
        help="Example: chr19 or comma-separated list chr17,chr18,chr19 or BED file or VCF file")

    main.add_argument("--threads", type=int, default=4,
        help="Threads for HaplotypeCaller and GenomicsDBImport (default: 4)")

    # -------------------------
    # Advanced arguments
    # -------------------------
    adv = parser.add_argument_group("Advanced arguments")

    adv.add_argument("--min-dp", type=int, default=20,
        help="Minimum read depth per SNP per sample (default: 20)")

    adv.add_argument("--min-qd", type=float, default=2.0,
        help="QD filter cutoff (default: 2.0)")

    adv.add_argument("--min-qual", type=float, default=30.0,
        help="QUAL filter cutoff (default: 30)")

    adv.add_argument("--java_mem", default="-Xmx32G",
        help="Java memory setting passed to GATK (default: -Xmx32G)")

    adv.add_argument("--tmp_dir", default=None,
        help="Directory for GATK temporary files (default: <output_dir>/tmp)")

    adv.add_argument("--platform", default="ILLUMINA",
        help="Read group platform (RGPL) used in AddOrReplaceReadGroups (default: ILLUMINA)")

    adv.add_argument("--call_SplitNCigar_for_assay_types",
        default="RNA,RNA-seq,RNAseq",
        help="Comma-separated assay names (case-insensitive) for which SplitNCigarReads will be run (default: RNA,RNA-seq,RNAseq)")

    adv.add_argument("--skip_AddOrReplaceRG", action="store_true",
        help="Skip AddOrReplaceReadGroups step (default: False). Caution: Use only if BAMs already contain correct RG fields.")

    adv.add_argument("--skip_SplitNCigar", action="store_true",
        help="Skip SplitNCigarReads for RNA-seq / spliced data (default: False).")

    adv.add_argument("--gatk_path", default="gatk",
        help="Path to GATK executable (default: gatk). Example: --gatk_path /path/to/gatk")

    adv.add_argument("--bcftools_path", default="bcftools",
        help="Path to bcftools executable (default: bcftools). Example: --bcftools_path /path/to/bcftools")

    adv.add_argument("--samtools_path", default="samtools",
        help="Path to samtools executable (default: samtools). Example: --samtools_path /path/to/samtools")

    # -------------------------
    # Slurm options for HaplotypeCaller
    # -------------------------
    adv.add_argument("--slurm_haplotypecaller", action="store_true",
        help="Run per-sample HaplotypeCaller on Slurm using sbatch (default: False).")

    adv.add_argument("--slurm_threads", type=int, default=None,
        help="Threads requested from Slurm for each HaplotypeCaller job (default: same as --threads).")

    adv.add_argument("--slurm_ram", default="32G",
        help="RAM requested from Slurm for each HaplotypeCaller job (default: 32G). Example: 16G or 64000M")

    adv.add_argument(
        "--slurm_pairwise_bcftools",
        action="store_true",
        help="Run pairwise bcftools queries on Slurm (default: False)"
    )
    
    # Internal-only worker flags (do not use directly)
    adv.add_argument("--_hc_slurm_worker", action="store_true", help=argparse.SUPPRESS)
    adv.add_argument("--_hc_sample_id", default=None, help=argparse.SUPPRESS)
    adv.add_argument("--_hc_processed_manifest", default=None, help=argparse.SUPPRESS)

    adv.add_argument("--_bcf_slurm_worker", action="store_true", help=argparse.SUPPRESS)
    adv.add_argument("--_bcf_pair", default=None, help=argparse.SUPPRESS)
    adv.add_argument("--_bcf_out", default=None, help=argparse.SUPPRESS)
    adv.add_argument("--_bcf_pass_vcf", default=None, help=argparse.SUPPRESS)
    return parser.parse_args()


def build_hc_slurm_script(sid, slurm_threads, slurm_ram, log_path, worker_cmd):
    return f"""#!/bin/bash
#SBATCH --job-name=HC_{sid}
#SBATCH --cpus-per-task={slurm_threads}
#SBATCH --mem={slurm_ram}
#SBATCH --output={log_path}

set -euo pipefail

echo "[info] Running HaplotypeCaller for {sid}"
echo "[info] Host: $(hostname)"
echo "[info] Start: $(date)"

{worker_cmd}

echo "[info] Done: $(date)"
"""


def build_bcf_slurm_script(s1, s2, log_path, worker_cmd):
    return f"""#!/bin/bash
#SBATCH --job-name=BCF_{s1}_{s2}
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output={log_path}

set -euo pipefail

{worker_cmd}
"""


def submit_and_wait_hc(hc_jobs, poll_seconds=30):
    job_ids = []
    for sid, script_path in hc_jobs:
        job_id = slurm_submit(script_path)
        print(f"[info] Submitted HaplotypeCaller job for {sid}: {job_id}")
        job_ids.append(job_id)
    slurm_wait(job_ids, poll_seconds=poll_seconds)


def submit_and_wait(job_scripts, poll_seconds=30):
    job_ids = []
    for script_path in job_scripts:
        job_id = slurm_submit(script_path)
        job_ids.append(job_id)
    slurm_wait(job_ids, poll_seconds=poll_seconds)



def main():
    ###############################
    # Argument parsing
    ###############################
    args = parse_args()
    
    ###############################
    # Worker dispatch
    ###############################
    if args._hc_slurm_worker:
        if args._hc_sample_id is None or args._hc_processed_manifest is None:
            raise ValueError("HC Slurm worker requires sample id and manifest")
        args.slurm_sample_id = args._hc_sample_id
        args.processed_manifest = args._hc_processed_manifest
        hc_slurm_worker(args)
        return
    
    if args._bcf_slurm_worker:
        if args._bcf_pair is None or args._bcf_out is None or args._bcf_pass_vcf is None:
            raise ValueError(
                "bcftools Slurm worker requires --_bcf_pair, --_bcf_out, and --_bcf_pass_vcf"
            )

        args._pair = args._bcf_pair
        args._pair_out = args._bcf_out
        args.pass_vcf = args._bcf_pass_vcf

        pairwise_bcftools_slurm_worker(args)
        return

    ###############################
    # Tool paths and timing
    ###############################
    tools = {
        "gatk": args.gatk_path,
        "bcftools": args.bcftools_path,
        "samtools": args.samtools_path,
    }

    # start timing
    start_time = time.time()
    start_timestamp = datetime.datetime.now()

    # normalize list of assays for SplitNCigarReads
    split_assays = [
        x.strip().lower()
        for x in args.call_SplitNCigar_for_assay_types.split(",")
        if x.strip()
    ]

    ###############################
    # Input loading and validation
    ###############################
    df = load_and_validate_input(args, tools["samtools"])

    ###############################
    # Output directory setup
    ###############################
    (
        outdir,
        tmp_dir,
        gvcf_dir,
        db_dir,
        joint_dir,
        rg_bam_dir,
        split_bam_dir,
    ) = prepare_output_dirs(args)
    
    gatk_java_opts = f'{args.java_mem} -Djava.io.tmpdir={tmp_dir}'
    
    ########################
    # Provenance log (start)
    ########################
    write_run_info_start(outdir, args, split_assays, start_timestamp)
    
    ###############################
    # Region arguments
    ###############################
    regions_arg = build_regions_arg(args.regions)

    ###############################
    # Preprocessing (RG + SplitNCigar)
    ###############################
    df, processed_manifest = process_bams(
        df=df,
        args=args,
        tools=tools,
        gatk_java_opts=gatk_java_opts,
        split_assays=split_assays,
        rg_bam_dir=rg_bam_dir,
        split_bam_dir=split_bam_dir,
        joint_dir=joint_dir,
        regions_arg=regions_arg,
    )

    ###############################
    # GVCF calling
    ###############################
    df["gvcf"] = None

    # Slurm settings
    slurm_threads = args.slurm_threads if args.slurm_threads is not None else args.threads

    if args.slurm_haplotypecaller:
        slurm_dir = outdir / "slurm_haplotypecaller"
        slurm_scripts = slurm_dir / "scripts"
        slurm_logs = slurm_dir / "logs"
        slurm_scripts.mkdir(exist_ok=True, parents=True)
        slurm_logs.mkdir(exist_ok=True, parents=True)

        # job_ids = []
        hc_jobs = []

        script_self = Path(__file__).resolve()

        for _, row in df.iterrows():
            sid = row["sample_id"]
            gvcf = gvcf_dir / f"{sid}.g.vcf.gz"
            df.loc[df.sample_id == sid, "gvcf"] = str(gvcf)

            if gvcf.exists():
                print(f"[warning] GVCF already exists for {sid}: {gvcf}")
                print("          Skipping HaplotypeCaller for this sample.")
                continue

            script_path = slurm_scripts / f"hc_{sid}.sbatch"
            log_path = slurm_logs / f"hc_{sid}.%j.out"

            # Build worker call: reuse same args as the driver run
            # We only add the internal worker selectors.
            cmdline = " ".join([subprocess.list2cmdline([x]) for x in os.sys.argv[1:]])

            worker_cmd = (
                f"python {script_self} {cmdline} "
                f"--_hc_slurm_worker "
                f"--_hc_sample_id {sid} "
                f"--_hc_processed_manifest {processed_manifest}"
            )

            script_text = build_hc_slurm_script(
                sid=sid,
                slurm_threads=slurm_threads,
                slurm_ram=args.slurm_ram,
                log_path=log_path,
                worker_cmd=worker_cmd,
            )

            with open(script_path, "w") as f:
                f.write(script_text)

            # job_id = slurm_submit(script_path)
            # print(f"[info] Submitted HaplotypeCaller job for {sid}: {job_id}")
            # job_ids.append(job_id)
            hc_jobs.append((sid, script_path))

        # Wait for all submitted jobs
        submit_and_wait_hc(hc_jobs, poll_seconds=30)

    else:
        for _, row in df.iterrows():
            sid = row["sample_id"]
            bam = row["processed_bam"]
            gvcf = gvcf_dir / f"{sid}.g.vcf.gz"
            df.loc[df.sample_id == sid, "gvcf"] = str(gvcf)

            if gvcf.exists():
                print(f"[warning] GVCF already exists for {sid}: {gvcf}")
                print("          Skipping HaplotypeCaller for this sample.")
                continue

            cmd = f"""
{tools["gatk"]} --java-options "{gatk_java_opts}" HaplotypeCaller \
  -R {args.reference} \
  -I {bam} \
  --emit-ref-confidence GVCF \
  -O {gvcf} {regions_arg} \
  --native-pair-hmm-threads {args.threads}
"""
            run(cmd)

    ###############################
    # 2. Sample map
    ###############################
    samplemap = joint_dir / "samplemap.txt"
    with open(samplemap, "w") as f:
        for _, row in df.iterrows():
            f.write(f"{row['sample_id']}\t{row['gvcf']}\n")

    ###############################
    # 3. GenomicsDBImport
    ###############################
    if (db_dir / "callset.json").exists():
        print(f"[warning] GenomicsDB workspace already exists: {db_dir}")
        print("          Skipping GenomicsDBImport.")
    else:
        cmd = f"""
{tools["gatk"]} --java-options "{gatk_java_opts}" GenomicsDBImport \
  -R {args.reference} \
  --sample-name-map {samplemap} \
  --genomicsdb-workspace-path {db_dir} {regions_arg} \
  --reader-threads {args.threads}
"""
        run(cmd)

    ###############################
    # 4. Joint calling
    ###############################
    joint_vcf = joint_dir / "joint_calls.vcf.gz"

    if joint_vcf.exists():
        print(f"[warning] Joint VCF already exists, skipping GenotypeGVCFs: {joint_vcf}")
    else:
        cmd = f"""
{tools["gatk"]} --java-options "{gatk_java_opts}" GenotypeGVCFs \
  -R {args.reference} \
  -V gendb://{db_dir} \
  -O {joint_vcf}
"""
        run(cmd)

    ###############################
    # 5. SNP-only
    ###############################
    snp_vcf = joint_dir / "joint_snps.vcf.gz"

    if snp_vcf.exists():
        print(f"[warning] SNP VCF already exists, skipping SelectVariants: {snp_vcf}")
    else:
        cmd = f"""
{tools["gatk"]} --java-options "{gatk_java_opts}" SelectVariants \
  -R {args.reference} \
  -V {joint_vcf} \
  --select-type-to-include SNP \
  -O {snp_vcf}
"""
        run(cmd)

    ###############################
    # 6. Filter + PASS-only
    ###############################
    filt_vcf = joint_dir / "joint_snps_filtered.vcf.gz"
    pass_vcf = joint_dir / "joint_snps_pass.vcf.gz"
    pass_vcf_csi = pass_vcf.with_suffix(pass_vcf.suffix + ".csi")
    pass_vcf_tbi = pass_vcf.with_suffix(pass_vcf.suffix + ".tbi")

    if pass_vcf.exists() and (pass_vcf_csi.exists() or pass_vcf_tbi.exists()):
        print(f"[warning] PASS VCF already exists, reusing: {pass_vcf}")
    else:
        cmd = f"""
{tools["gatk"]} --java-options "{gatk_java_opts}" VariantFiltration \
  -R {args.reference} \
  -V {snp_vcf} \
  --filter-expression "QD < {args.min_qd} || QUAL < {args.min_qual}" \
  --filter-name LowQual \
  -O {filt_vcf}
{tools["bcftools"]} view -f PASS {filt_vcf} -Oz -o {pass_vcf}
{tools["bcftools"]} index {pass_vcf}
"""
        run(cmd)
    ###############################
    # Pairwise mismatch computation
    ###############################
    results = []

    if args.slurm_pairwise_bcftools:
        script_self = Path(__file__).resolve()
        cmdline = " ".join([subprocess.list2cmdline([x]) for x in os.sys.argv[1:]])

        pair_tmp_dir = joint_dir / PAIRWISE_TMP_DIR
        pair_tmp_dir.mkdir(exist_ok=True, parents=True)

        slurm_dir = outdir / "slurm_pairwise_bcftools"
        scripts_dir = slurm_dir / "scripts"
        logs_dir = slurm_dir / "logs"
        scripts_dir.mkdir(exist_ok=True, parents=True)
        logs_dir.mkdir(exist_ok=True, parents=True)

        # job_ids = []
        bcf_scripts = []

        for s1, s2 in combinations(df["sample_id"], 2):
            out_path = pair_tmp_dir / f"{s1}__{s2}.tsv"
            if out_path.exists():
                continue

            script_path = scripts_dir / f"bcf_{s1}__{s2}.sbatch"
            log_path = logs_dir / f"bcf_{s1}__{s2}.%j.out"

            worker_cmd = (
                f"python {script_self} {cmdline} "
                f"--_bcf_slurm_worker "
                f"--_bcf_pair {s1},{s2} "
                f"--_bcf_out {out_path} "
                f"--_bcf_pass_vcf {pass_vcf}"
            )

            script_text = build_bcf_slurm_script(
                s1=s1,
                s2=s2,
                log_path=log_path,
                worker_cmd=worker_cmd,
            )

            script_path.write_text(script_text)

            bcf_scripts.append(script_path)
        submit_and_wait(bcf_scripts, poll_seconds=30)
        
        for s1, s2 in combinations(df["sample_id"], 2):
            out_path = pair_tmp_dir / f"{s1}__{s2}.tsv"
            if not out_path.exists():
                continue

            lines = out_path.read_text().splitlines()

            total, mismatches, perc = parse_bcftools_pairwise_output(
                lines=lines,
                min_dp=args.min_dp,
                min_sites_warning=MIN_SITES_WARNING,
                s1=s1,
                s2=s2,
            )
            
            results.append([s1, s2, total, mismatches, perc])
    else:
        # ---------- Serial execution ----------
        for s1, s2 in combinations(df["sample_id"], 2):

            cmd = build_bcftools_pairwise_query(
                bcftools=tools["bcftools"],
                pass_vcf=pass_vcf,
                s1=s1,
                s2=s2,
            )

            lines = run_bcftools_pairwise_query(cmd)

            total, mismatches, perc = parse_bcftools_pairwise_output(
                lines=lines,
                min_dp=args.min_dp,
                min_sites_warning=MIN_SITES_WARNING,
                s1=s1,
                s2=s2,
            )

            results.append([s1, s2, total, mismatches, perc])
    
    res = pd.DataFrame(
        results,
        columns=["sample1", "sample2", "total_sites", "mismatches", "percent_mismatch"]
    )
    ###############################
    # 8. Nearest-neighbor alerts
    ###############################

    # map sample -> assay
    sample_to_assay = dict(zip(df["sample_id"], df["assay"]))

    ###############################
    # Nearest-neighbor reporting
    ###############################
    print("\n[Nearest-neighbor genotype similarity alerts]\n")

    different_rows = []
    same_rows = []
    
    # --- make pairwise table symmetric for NN only ---
    res_nn = pd.concat(
        [
            res,
            res.rename(columns={"sample1": "sample2", "sample2": "sample1"})
        ],
        ignore_index=True
    )
    
    # for s in res["sample1"].unique():
    for s in res_nn["sample1"].unique():

        assay_s = sample_to_assay.get(s, "NA")

        # sub = res[res["sample1"] == s].dropna(subset=["percent_mismatch"])
        sub = res_nn[res_nn["sample1"] == s].dropna(subset=["percent_mismatch"])
        if sub.empty:
            continue

        ##########################
        # Different-assay neighbor
        ##########################
        cross = sub[sub["sample2"].map(lambda x: sample_to_assay.get(x, None) != assay_s)]

        if not cross.empty:
            best_cross = cross.sort_values("percent_mismatch").head(1).iloc[0]

            print(
                f"[Different assay] Closest match for {s} "
                f"(assay={assay_s}) is {best_cross['sample2']} "
                f"(assay={sample_to_assay.get(best_cross['sample2'], 'NA')}) "
                f"({best_cross['percent_mismatch']:.2f}%)"
            )

            different_rows.append({
                "sample": s,
                "sample_assay": assay_s,
                "closest_sample": best_cross["sample2"],
                "closest_sample_assay": sample_to_assay.get(best_cross["sample2"], "NA"),
                "percent_mismatch": best_cross["percent_mismatch"]
            })

        ##########################
        # Same-assay neighbor
        ##########################
        same = sub[sub["sample2"].map(lambda x: sample_to_assay.get(x, None) == assay_s)]
        same = same[same["sample2"] != s]  # exclude self

        if not same.empty:
            best_same = same.sort_values("percent_mismatch").head(1).iloc[0]

            print(
                f"[Same assay] Closest match for {s} "
                f"(assay={assay_s}) is {best_same['sample2']} "
                f"({best_same['percent_mismatch']:.2f}%)"
            )

            same_rows.append({
                "sample": s,
                "sample_assay": assay_s,
                "closest_sample": best_same["sample2"],
                "closest_sample_assay": assay_s,
                "percent_mismatch": best_same["percent_mismatch"]
            })

    ###############################
    # Save nearest neighbor report
    ###############################
    nn_path = outdir / "nearest_neighbor_matches.tsv"

    with open(nn_path, "w") as f:

        f.write("### Closest match from a different assay\n")
        f.write("sample\tsample_assay\tclosest_sample\tclosest_sample_assay\tpercent_mismatch\n")
        for r in different_rows:
            f.write(
                f"{r['sample']}\t{r['sample_assay']}\t"
                f"{r['closest_sample']}\t{r['closest_sample_assay']}\t"
                f"{r['percent_mismatch']:.6f}\n"
            )

        f.write("\n### Closest match from the same assay\n")
        f.write("sample\tsample_assay\tclosest_sample\tclosest_sample_assay\tpercent_mismatch\n")
        for r in same_rows:
            f.write(
                f"{r['sample']}\t{r['sample_assay']}\t"
                f"{r['closest_sample']}\t{r['closest_sample_assay']}\t"
                f"{r['percent_mismatch']:.6f}\n"
            )

    print(f"\nWrote nearest-neighbor summary: {nn_path}\n")

    ###############################
    # Final outputs and provenance
    ###############################
    tsv_path = outdir / "genotype_mismatch_rates.tsv"
    csv_path = outdir / "genotype_mismatch_rates.csv"

    res["percent_mismatch"] = res["percent_mismatch"].round(5)

    res.to_csv(tsv_path, sep="\t", index=False)
    res.to_csv(csv_path, index=False)

    
    write_run_info_end(outdir, start_time)
    
    print("\nDone.")
    print(f"Wrote:\n  {tsv_path}\n  {csv_path}\n")


if __name__ == "__main__":
    main()
