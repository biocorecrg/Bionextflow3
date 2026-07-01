#!/usr/bin/env python

import pandas as pd
import argparse
import sys

# Usage: python smallRNA_summary.py desc.txt annotated_counts.txt
parser = argparse.ArgumentParser(description="Creates table for each type of small RNA counts.")
parser.add_argument("-d","--desc", help="Path to the desc.txt file.")
parser.add_argument("-a","--annotation", help="Path to the annotated_counts.txt file.")
parser.add_argument("-e","--experiment", help='Indicate experiment data comes from: "rnaseq" or "smallrnaseq"')
parser.add_argument("-r","--biotype", type=str, help="Comma-separated list of RNA types to display in the bargraph (e.g. 'miRNA,snRNA')")
args = parser.parse_args()


desc_file = args.desc
annotated_file = args.annotation
experiment = args.experiment

# Set small_RNAs default based on experiment if not provided
if args.biotype:
    biotypes = [rna.strip() for rna in args.biotype.split(",") if rna.strip()]
elif experiment == "rnaseq":
    biotypes = ["protein_coding", "lncRNA","rRNA"]
else:
    biotypes = ["miRNA","snoRNA","rRNA","protein_coding","intergenic"] 

# --- check files exist ---
for f in [desc_file, annotated_file]:
    try:
        open(f).close()
    except FileNotFoundError:
        print(f"File {f} does not exist.")
        sys.exit(1)

# --- load description file (sample names are in first column, skip header) ---
desc_df = pd.read_csv(desc_file, sep="\t", usecols=[0])
sample_names = desc_df.iloc[0:, 0].tolist()  # skip header

# --- load annotated counts file ---

# Load annotated counts file

anno_df = pd.read_csv(annotated_file, sep=",")

# Normalize gene.type column:
# 1. intron-* types are collapsed to "intron"
# 2. When mixing feature + intron, the feature always wins
# 3. For smallRNA experiments: smallRNA types beat protein_coding
# 4. Deduplicate remaining types
SMALL_RNA_TYPES = {"miRNA", "snoRNA", "snRNA", "piRNA"}

def normalize_gene_type(val):
    if pd.isna(val):
        return val
    parts = [p.strip() for p in str(val).split(",")]

    # Collapse intron-* to "intron"
    normalized = ["intron" if p.startswith("intron-") else p for p in parts]

    # Deduplicate while preserving order
    unique = list(dict.fromkeys(normalized))

    # If mix of intron and non-intron features, drop intron (feature wins)
    non_intron = [u for u in unique if u != "intron"]
    if non_intron and "intron" in unique:
        unique = non_intron

    # For smallRNA experiments: smallRNA types take priority over protein_coding
    if experiment != "rnaseq" and "protein_coding" in unique:
        small_rnas = [u for u in unique if u in SMALL_RNA_TYPES]
        if small_rnas:
            unique = [u for u in unique if u != "protein_coding"]

    if len(unique) == 1:
        return unique[0]

    # Still ambiguous (multiple different features)
    return ",".join(unique)

anno_df["gene.type"] = anno_df["gene.type"].apply(normalize_gene_type)


rows = []

for sample in sample_names:
    if sample not in anno_df.columns:
        print(f"Warning: sample {sample} not found in annotated file.")
        continue

    ## selecting the gene.type column and the column with the counts for the selected sample
    subsample = anno_df[["gene.type", sample]].copy()
    subsample.columns = ["gene.type", "count"]

    # ambiguous counts (gene.type contains ',' and all values are different) 
    ambiguous_counts = subsample.loc[subsample["gene.type"].str.contains(",", na=False), "count"].sum()

    subsample.loc[subsample["gene.type"].str.contains(",", na=False), "gene.type"]

    # non-ambiguous counts (gene.type does not contain ',')
    non_ambiguous = subsample.loc[~subsample["gene.type"].str.contains(",", na=False), "count"]
    non_ambiguous_counts = non_ambiguous.sum()

    # total counts
    total_counts = subsample["count"].sum()

    # counts for each RNA type
    counts = []
    sum_small = 0
    for small in biotypes:
        small_count = subsample.loc[(subsample["gene.type"] == small) & ~subsample["gene.type"].str.contains(",", na=False), "count"].sum()
        counts.append(small_count)
        sum_small += small_count

    other = non_ambiguous_counts - sum_small

    row = [sample] + counts + [ambiguous_counts, other]
    rows.append(row)


# --- save output ---
columns = ["sample_name"] + biotypes + ["ambiguous", "other"]
summary_df = pd.DataFrame(rows, columns=columns)
summary_df.to_csv("RNAstats_mqc.csv", sep=",", index=False)
