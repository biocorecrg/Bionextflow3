#!/usr/bin/env python3
import sys
import os
import glob
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Extract the accepted reads from an adaptive sampling experiment."
    )
    parser.add_argument(
        "-a", "--as_decisions",
        nargs="+",
        required=True,
        help="Path to one or more adaptive sampling decisions files."
    )
    parser.add_argument(
        "-s", "--sequencing_summary",
        nargs="+",
        required=True,
        help="Path to one or more sequencing summary files."
    )
    parser.add_argument(
        "-o", "--output_name",
        required=True,
        help="Base name for the output file (generates <output_name>-Accepted-Reads.txt)."
    )
    
    args = parser.parse_args()
    output_name = args.output_name

    # 1. Read accepted read IDs from adaptive sampling decisions
    sequenced_reads = set()
    for file_path in args.as_decisions:
        try:
            with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    parts = line.strip().split(",")
                    if len(parts) >= 2 and parts[1] == "sequence":
                        sequenced_reads.add(parts[0])
        except Exception as e:
            print(f"Warning: Could not read {file_path}: {e}", file=sys.stderr)

    if not sequenced_reads:
        print("No sequencing reads identified in adaptive sampling decisions.", file=sys.stderr)
    
    # 2. Filter read IDs from sequencing summary files
    output_file_path = f"{output_name}-Accepted-Reads.txt"
    
    total_reads = 0
    accepted_reads = 0

    try:
        with open(output_file_path, "w", encoding="utf-8") as out_f:
            for file_path in args.sequencing_summary:
                try:
                    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
                        is_first_line = True
                        for line in f:
                            parts = line.strip().split("\t")
                            if len(parts) >= 5:
                                if is_first_line and parts[4] == "parent_read_id":
                                    is_first_line = False
                                    continue
                                read_id = parts[4]
                                total_reads += 1
                                if read_id in sequenced_reads:
                                    accepted_reads += 1
                                    out_f.write(read_id + "\n")
                except Exception as e:
                    print(f"Warning: Could not read {file_path}: {e}", file=sys.stderr)
    except Exception as e:
        print(f"Error: Could not write output to {output_file_path}: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Successfully wrote accepted reads to {output_file_path}")

    percentage = (accepted_reads / total_reads * 100) if total_reads > 0 else 0.0
    stats_file_path = f"{output_name}_mqc.txt"
    try:
        with open(stats_file_path, "w", encoding="utf-8") as stat_f:
            stat_f.write("# format: 'tsv'\n")
            stat_f.write("# plot_type: 'table'\n")
            stat_f.write("# section_name: 'Adaptive Sampling Reads'\n")
            stat_f.write("# id: 'adaptive_sampling_reads'\n")
            stat_f.write("# description: 'Statistics of reads from the adaptive sampling decisions.'\n")
            stat_f.write("# pconfig:\n")
            stat_f.write("#    namespace: 'Adaptive Sampling'\n")
            stat_f.write(f"Sample\tTotal Reads\tAccepted Reads\t% Accepted\n")
            stat_f.write(f"{output_name}\t{total_reads}\t{accepted_reads}\t{percentage:.2f}\n")
        print(f"Successfully wrote MultiQC stats to {stats_file_path}")
    except Exception as e:
        print(f"Warning: Could not write stats file {stats_file_path}: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()
