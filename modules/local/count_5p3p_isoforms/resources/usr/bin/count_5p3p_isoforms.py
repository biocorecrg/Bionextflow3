#!/usr/bin/env python
"""
Count reads per known miRNA isoform (5p/3p) from ShortStack condensed BAMs.

Uses known_miRNAs.gff3 from ShortStack for miRNA coordinates and names
(which already include -5p/-3p suffixes from miRBase).
Reads the XW:i tag from condensed BAMs for true read abundance.
Reads must not be larger than the miRNA length plus one base.

Outputs a single matrix with miRNA, followed by 5p and 3p fractions
(A/(A+B) and B/(A+B)) for each sample.
"""

import pysam
import csv
import argparse
import os
import re
from collections import defaultdict


def parse_gff3(gff3_file):
    """Parse known_miRNAs.gff3 to get miRNA regions.

    Returns:
        regions: list of dict of {chrom, start, end, strand, name}
    """
    coord_to_names = defaultdict(list)
    seen_coords = []

    with open(gff3_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom = fields[0]
            start = int(fields[3])  # 1-based start
            end = int(fields[4])    # 1-based end
            strand = fields[6]

            # Parse Name from attributes
            name = None
            for attr in fields[8].split(';'):
                attr = attr.strip()
                if attr.startswith('Name='):
                    name = attr.split('=', 1)[1]

            if name:
                coord = (chrom, start, end, strand)
                if coord not in coord_to_names:
                    seen_coords.append(coord)
                if name not in coord_to_names[coord]:
                    coord_to_names[coord].append(name)

    regions = []
    for coord in seen_coords:
        chrom, start, end, strand = coord
        names_list = coord_to_names[coord]
        merged_name = ",".join(names_list)
        regions.append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'strand': strand,
            'name': merged_name
        })

    return regions


def count_reads_in_bam(bam_file, regions, strandness='stranded'):
    """Count reads from a condensed BAM overlapping each miRNA region.

    Returns dict: (chrom, start, end, name) -> total count
    """
    # Index regions by (chrom, strand) for fast lookup
    regions_by_key = defaultdict(list)
    for r in regions:
        key = (r['chrom'], r['strand'])
        # Store 0-based coordinates for pysam overlap matching
        regions_by_key[key].append({
            'start_0': r['start'] - 1,
            'end_0': r['end'],
            'orig_r': r
        })

    counts = defaultdict(int)

    bam = pysam.AlignmentFile(bam_file, 'rb')

    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            continue

        read_chrom = read.reference_name
        bam_strand = '-' if read.is_reverse else '+'
        
        if strandness == 'stranded':
            read_strands = [bam_strand]
        elif strandness == 'reverse':
            read_strands = ['+' if bam_strand == '-' else '-']
        elif strandness == 'unstranded':
            read_strands = ['+', '-']
        else:
            read_strands = [bam_strand]

        read_start = read.reference_start    # 0-based
        read_end = read.reference_end         # 0-based, exclusive

        for r_strand in read_strands:
            key = (read_chrom, r_strand)
            if key not in regions_by_key:
                continue

            # Get read count from XW tag (condensed read depth)
            try:
                xw = read.get_tag('XW')
            except KeyError:
                xw = 1  # Fallback if tag missing

            # Check overlap with miRNA regions on same chrom/strand
            for r_info in regions_by_key[key]:
                if read_start < r_info['end_0'] and read_end > r_info['start_0']:
                    region_len = r_info['end_0'] - r_info['start_0']
                    read_len = read.query_length if read.query_length else (read_end - read_start)
                    if read_len <= region_len + 1:
                        orig = r_info['orig_r']
                        r_key = (orig['chrom'], orig['start'], orig['end'], orig['name'])
                        counts[r_key] += xw

    bam.close()
    return counts


def get_sample_name(bam_path):
    """Extract clean sample name from BAM filename."""
    name = os.path.basename(bam_path)
    # Remove common suffixes added by the pipeline
    name = re.sub(r'_trimmed_condensed\.bam$', '', name)
    name = re.sub(r'_filter\.bam$', '', name)
    name = re.sub(r'\.bam$', '', name)
    return name


def main():
    parser = argparse.ArgumentParser(
        description="Count reads per known miRNA isoform from ShortStack condensed BAMs"
    )
    parser.add_argument(
        "--gff3", required=True,
        help="known_miRNAs.gff3 from ShortStack"
    )
    parser.add_argument(
        "--bams", required=True, nargs='+',
        help="Per-sample condensed BAM files"
    )
    parser.add_argument(
        "--output", required=True,
        help="Output TSV fraction matrix"
    )
    parser.add_argument(
        "--output-counts",
        help="Optional output TSV file for raw read counts"
    )
    parser.add_argument(
        "--strandness", default="stranded",
        choices=["stranded", "reverse", "unstranded"],
        help="Strand specificity mode for matching reads to miRNAs"
    )

    args = parser.parse_args()

    # Parse miRNA regions from GFF3
    regions = parse_gff3(args.gff3)

    if not regions:
        # No miRNAs found, write empty output with header only
        with open(args.output, 'w') as f:
            f.write("miRNA\n")
        if args.output_counts:
            with open(args.output_counts, 'w') as f:
                f.write("miRNA\n")
        return

    # Identify paired miRNAs globally
    base_to_arms = defaultdict(set)
    for r in regions:
        for name in r['name'].split(','):
            if name.endswith('-5p'):
                base_to_arms[name[:-3]].add('5p')
            elif name.endswith('-3p'):
                base_to_arms[name[:-3]].add('3p')

    paired_miRNAs = sorted([base for base, arms in base_to_arms.items() if '5p' in arms and '3p' in arms])

    # Count reads per sample
    sample_names = []
    sample_counts = {}

    for bam_file in sorted(args.bams):
        sample_name = get_sample_name(bam_file)
        sample_names.append(sample_name)
        sample_counts[sample_name] = count_reads_in_bam(bam_file, regions, args.strandness)

    # Write output matrix
    # Format: miRNA, sample1_5p, sample1_3p, sample2_5p, sample2_3p, ...
    with open(args.output, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        header = ['miRNA']
        for sample in sample_names:
            header.extend([f"{sample}_5p", f"{sample}_3p"])
        writer.writerow(header)

        for base in paired_miRNAs:
            row = [base]
            for sample in sample_names:
                # Sum all coordinates for base-5p globally
                sum_5p = sum(
                    val for r_key, val in sample_counts[sample].items()
                    if f"{base}-5p" in r_key[3].split(',')
                )
                # Sum all coordinates for base-3p globally
                sum_3p = sum(
                    val for r_key, val in sample_counts[sample].items()
                    if f"{base}-3p" in r_key[3].split(',')
                )

                total = sum_5p + sum_3p
                if total > 0:
                    prop_5p = round(sum_5p / total, 4)
                    prop_3p = round(sum_3p / total, 4)
                else:
                    prop_5p = 0.0
                    prop_3p = 0.0

                row.extend([prop_5p, prop_3p])
            writer.writerow(row)

    # Write raw counts matrix if requested
    if args.output_counts:
        all_miRNAs = sorted(list(set(r['name'] for r in regions)))
        with open(args.output_counts, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            header = ['miRNA'] + sample_names
            writer.writerow(header)

            for name in all_miRNAs:
                # Format row name: replace -5p / -3p with _5p / _3p
                row_name = name.replace('-5p', '_5p').replace('-3p', '_3p')
                row = [row_name]
                for sample in sample_names:
                    total_count = sum(
                        val for r_key, val in sample_counts[sample].items()
                        if r_key[3] == name
                    )
                    row.append(total_count)
                writer.writerow(row)


if __name__ == "__main__":
    main()
