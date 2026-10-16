#!/usr/bin/env python3
"""
Wrapper for bam2baseQ.py:
1. Pre-checks input BAM/CRAM files using pysam.idxstats. If mapped read count is below min_cov,
   generates an empty BED (.bed or .bed.gz) and exits cleanly.
2. Runs /opt/app/src/bam2baseQ.py.
3. Catches crashes due to zero covered regions ("0 bases in 0 regions to process" /
   "ValueError: Number of processes must be at least 1") and cleanly outputs an empty BED file.
"""

import sys
import os
import gzip
import subprocess
import argparse
import pysam

def get_total_mapped_reads(bam_path):
    try:
        stats = pysam.idxstats(bam_path)
        total_mapped = sum(int(line.split('\t')[2]) for line in stats.splitlines() if line.strip())
        return total_mapped
    except Exception as e:
        sys.stderr.write(f"[run_bam2baseQ] Warning: could not retrieve idxstats for {bam_path}: {e}\n")
        return None

def write_empty_output(out_path):
    if not out_path or out_path == "-":
        return
    if out_path.endswith(".gz"):
        with gzip.open(out_path, "wb") as f:
            pass
    else:
        with open(out_path, "w") as f:
            pass

def main():
    parser = argparse.ArgumentParser(
        description="Wrapper for bam2baseQ.py with zero-coverage pre-check and error handling",
        add_help=False
    )
    parser.add_argument("-i", "--input", nargs="+", required=False, default=[])
    parser.add_argument("-o", "--output", default=None)
    parser.add_argument("-c", "--cov", type=int, default=10)
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("-h", "--help", action="store_true")
    parser.add_argument("--version", action="store_true")

    args, unknown_args = parser.parse_known_args()

    if args.help:
        cmd = ["/opt/app/src/bam2baseQ.py", "--help"]
        subprocess.run(cmd)
        sys.exit(0)

    if args.version:
        cmd = ["/opt/app/src/bam2baseQ.py", "--version"]
        subprocess.run(cmd)
        sys.exit(0)

    bams = args.input
    if len(bams) < 2:
        sys.stderr.write("[run_bam2baseQ] Error: At least two input BAM/CRAM files are required (-i file1 file2).\n")
        sys.exit(1)

    min_cov = args.cov
    out_file = args.output

    # 1. Pre-check mapped reads
    for b in bams:
        mapped = get_total_mapped_reads(b)
        if mapped is not None and mapped < min_cov:
            sys.stderr.write(
                f"[run_bam2baseQ] Info: Input alignment '{b}' has {mapped} mapped reads, "
                f"which is less than min coverage threshold ({min_cov}).\n"
                f"[run_bam2baseQ] Info: Skipping bam2baseQ and generating empty output.\n"
            )
            write_empty_output(out_file)
            sys.exit(0)

    # 2. Run /opt/app/src/bam2baseQ.py
    cmd = ["/opt/app/src/bam2baseQ.py", "-i", bams[0], bams[1], "-c", str(min_cov), "-t", str(args.threads)] + unknown_args

    is_gz_output = out_file and out_file.endswith(".gz")

    if is_gz_output:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        with gzip.open(out_file, "wt") as gz_out:
            for line in proc.stdout:
                gz_out.write(line)
        _, stderr_text = proc.communicate()
    elif out_file and out_file != "-":
        with open(out_file, "w") as f_out:
            proc = subprocess.run(cmd, stdout=f_out, stderr=subprocess.PIPE, text=True)
            stderr_text = proc.stderr
    else:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        for line in proc.stdout:
            sys.stdout.write(line)
        _, stderr_text = proc.communicate()

    if proc.returncode != 0:
        if "0 bases in 0 regions to process" in stderr_text or "ValueError: Number of processes must be at least 1" in stderr_text:
            sys.stderr.write(
                f"[run_bam2baseQ] Warning: Zero regions covered with >= {min_cov} reads. "
                f"Generating empty output file.\n"
            )
            write_empty_output(out_file)
            sys.exit(0)
        else:
            sys.stderr.write(stderr_text)
            sys.exit(proc.returncode)
    else:
        if stderr_text:
            sys.stderr.write(stderr_text)

if __name__ == "__main__":
    main()
