#!/usr/bin/env python3

import pandas as pd
import scipy.sparse as sp
import scipy.io as sio
import gzip
import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert STAR ReadsPerGene to Cell Ranger sparse matrix format"
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Directory with ReadsPerGene.out.tab files"
    )
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument(
        "-d",
        "--desc",
        required=True,
        help="Gene description file (gene_id, gene_name, gene_type)",
    )
    parser.add_argument(
        "-s",
        "--strand",
        default=2,
        type=int,
        choices=[1, 2, 3],
        help="Strand column: 1=unstranded, 2=sense, 3=antisense (default: 2)",
    )
    parser.add_argument(
        "-e",
        "--ext",
        default="ReadsPerGene.out.tab",
        help="File extension to look for (default: ReadsPerGene.out.tab)",
    )
    parser.add_argument(
        "--strip-version",
        action="store_true",
        default=True,
        help="Strip Ensembl version from gene IDs (default: True)",
    )
    return parser.parse_args()


def load_gene_map(desc_file):
    desc = pd.read_csv(desc_file, sep="\t", index_col=0)
    desc.index = desc.index.str.strip('"')
    desc["gene_name"] = desc["gene_name"].str.strip('"')
    return desc["gene_name"].to_dict()


def load_counts(input_dir, ext, strand_col):
    rows = []
    cols = []
    data = []
    cell_ids = []
    gene_ids = None
    gene_to_idx = {}

    files = sorted([f for f in os.listdir(input_dir) if f.endswith(ext)])
    if not files:
        raise ValueError(f"No files found with extension '{ext}' in {input_dir}")

    for col_idx, f in enumerate(files):
        cell_id = f.replace(ext, "").rstrip(".")
        cell_ids.append(cell_id)
        fpath = os.path.join(input_dir, f)
        
        # Read only the necessary column to save memory
        df = pd.read_csv(
            fpath, sep="\t", skiprows=4, header=None, usecols=[0, strand_col], index_col=0
        )
        
        if gene_ids is None:
            gene_ids = df.index.tolist()
            gene_to_idx = {gene: i for i, gene in enumerate(gene_ids)}
        
        series = df[strand_col]
        nonzero = series[series > 0]
        
        if not nonzero.empty:
            data.extend(nonzero.values)
            cols.extend([col_idx] * len(nonzero))
            rows.extend([gene_to_idx[g] for g in nonzero.index])

    matrix = sp.coo_matrix(
        (data, (rows, cols)), shape=(len(gene_ids), len(cell_ids)), dtype="int32"
    )
    return matrix.tocsr(), gene_ids, cell_ids


def write_output(matrix, gene_ids, cell_ids, gene_map, output_dir, strip_version):
    os.makedirs(output_dir, exist_ok=True)

    # barcodes.tsv.gz
    with gzip.open(os.path.join(output_dir, "barcodes.tsv.gz"), "wt") as f:
        for cell in cell_ids:
            f.write(cell + "\n")

    # features.tsv.gz
    with gzip.open(os.path.join(output_dir, "features.tsv.gz"), "wt") as f:
        for gene in gene_ids:
            gene_id = gene.split(".")[0] if strip_version else gene
            gene_name = gene_map.get(gene, gene_id)
            f.write(f"{gene_id}\t{gene_name}\tGene Expression\n")

    # matrix.mtx.gz
    with gzip.open(os.path.join(output_dir, "matrix.mtx.gz"), "wb") as f:
        sio.mmwrite(f, matrix)


def main():
    args = parse_args()

    print(f"Loading gene map from {args.desc}...")
    gene_map = load_gene_map(args.desc)

    print(f"Loading counts from {args.input}...")
    matrix, gene_ids, cell_ids = load_counts(args.input, args.ext, args.strand)

    print(f"Writing output to {args.output}...")
    write_output(matrix, gene_ids, cell_ids, gene_map, args.output, args.strip_version)

    print(f"Done: {len(gene_ids)} genes x {len(cell_ids)} cells")


if __name__ == "__main__":
    main()
