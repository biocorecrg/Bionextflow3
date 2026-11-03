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
    accepted_actions = {"sequence", "stop_receiving"}

    for file_path in args.as_decisions:
        try:
            with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    parts = line.strip().split(",")
                    if len(parts) >= 2 and parts[1].strip() in accepted_actions:
                        sequenced_reads.add(parts[0].strip())
        except Exception as e:
            print(f"Warning: Could not read {file_path}: {e}", file=sys.stderr)

    if not sequenced_reads:
        print("No sequencing reads identified in adaptive sampling decisions.", file=sys.stderr)

    # 2. Filter read IDs from sequencing summary files
    output_file_path = f"{output_name}-Accepted-Reads.txt"

    total_reads = 0
    accepted_reads = 0
    written_reads = set()

    try:
        with open(output_file_path, "w", encoding="utf-8") as out_f:
            for file_path in args.sequencing_summary:
                try:
                    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
                        read_id_idx = 3
                        parent_id_idx = 4

                        for line in f:
                            parts = line.strip().split("\t")
                            if not parts or len(parts) < 2:
                                continue

                            # Detect header row and map column indices dynamically
                            if "read_id" in parts:
                                read_id_idx = parts.index("read_id")
                                parent_id_idx = (
                                    parts.index("parent_read_id")
                                    if "parent_read_id" in parts
                                    else read_id_idx
                                )
                                continue

                            if len(parts) <= max(read_id_idx, parent_id_idx):
                                continue

                            read_id = parts[read_id_idx].strip()
                            parent_id = parts[parent_id_idx].strip()

                            total_reads += 1

                            # Match against AS decision using parent_read_id (or read_id),
                            # but always write the actual read_id present in the BAM/CRAM
                            if parent_id in sequenced_reads or read_id in sequenced_reads:
                                accepted_reads += 1
                                if read_id not in written_reads:
                                    out_f.write(read_id + "\n")
                                    written_reads.add(read_id)
                except Exception as e:
                    print(f"Warning: Could not read {file_path}: {e}", file=sys.stderr)
    except Exception as e:
        print(f"Error: Could not write output to {output_file_path}: {e}", file=sys.stderr)
        sys.exit(1)

    print(
        f"Successfully wrote {len(written_reads)} unique accepted reads to {output_file_path}"
    )

    # 3. Write MultiQC stats
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
            stat_f.write("Sample\tTotal Reads\tAccepted Reads\t% Accepted\n")
            stat_f.write(f"{output_name}\t{total_reads}\t{accepted_reads}\t{percentage:.2f}\n")
        print(f"Successfully wrote MultiQC stats to {stats_file_path}")
    except Exception as e:
        print(f"Warning: Could not write stats file {stats_file_path}: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
