#!/usr/bin/env python3

import numpy as np
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
    """Load gene_id -> gene_name mapping without pandas, line by line."""
    gene_map = {}
    with open(desc_file, "r") as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                gene_id = parts[0].strip('"')
                gene_name = parts[1].strip('"')
                gene_map[gene_id] = gene_name
    return gene_map


def discover_genes(input_dir, ext, strand_col):
    """Read one file to discover the gene order (after skipping the 4 header rows)."""
    for fname in sorted(os.listdir(input_dir)):
        if fname.endswith(ext):
            genes = []
            fpath = os.path.join(input_dir, fname)
            with open(fpath, "r") as f:
                for i, line in enumerate(f):
                    if i < 4:
                        continue
                    gene_id = line.split("\t", 1)[0]
                    genes.append(gene_id)
            return genes
    return []


def build_sparse_matrix(input_dir, ext, strand_col, gene_list):
    """
    Stream through files one at a time, collecting only non-zero entries
    into COO format arrays.  Never builds a dense matrix.
    """
    gene_to_idx = {g: i for i, g in enumerate(gene_list)}
    n_genes = len(gene_list)

    # COO accumulators
    rows = []
    cols = []
    vals = []

    barcodes = []
    col_idx = 0

    files = sorted(f for f in os.listdir(input_dir) if f.endswith(ext))
    if not files:
        raise ValueError(f"No files found with extension '{ext}' in {input_dir}")

    for fi, fname in enumerate(files):
        cell_id = fname.replace(ext, "").rstrip(".")
        barcodes.append(cell_id)

        fpath = os.path.join(input_dir, fname)
        with open(fpath, "r") as fh:
            for i, line in enumerate(fh):
                if i < 4:
                    continue
                parts = line.rstrip("\n").split("\t")
                gene_id = parts[0]
                count = int(parts[strand_col])
                if count != 0:
                    row = gene_to_idx.get(gene_id)
                    if row is not None:
                        rows.append(row)
                        cols.append(col_idx)
                        vals.append(count)

        col_idx += 1

        if (fi + 1) % 500 == 0:
            print(f"  Processed {fi + 1}/{len(files)} files...")

    # Build sparse matrix in one shot from COO data
    sparse_mtx = sp.coo_matrix(
        (np.array(vals, dtype=np.int32), (np.array(rows, dtype=np.int32), np.array(cols, dtype=np.int32))),
        shape=(n_genes, len(barcodes)),
    )

    return sparse_mtx, barcodes


def write_output(sparse_mtx, gene_list, gene_map, barcodes, output_dir, strip_version):
    os.makedirs(output_dir, exist_ok=True)

    # barcodes.tsv.gz
    with gzip.open(os.path.join(output_dir, "barcodes.tsv.gz"), "wt") as f:
        for cell in barcodes:
            f.write(cell + "\n")

    # features.tsv.gz
    with gzip.open(os.path.join(output_dir, "features.tsv.gz"), "wt") as f:
        for gene in gene_list:
            gene_id = gene.split(".")[0] if strip_version else gene
            gene_name = gene_map.get(gene, gene_id)
            f.write(f"{gene_id}\t{gene_name}\tGene Expression\n")

    # matrix.mtx.gz — write CSC (Cell Ranger convention)
    sparse_csc = sparse_mtx.tocsc()
    with gzip.open(os.path.join(output_dir, "matrix.mtx.gz"), "wb") as f:
        sio.mmwrite(f, sparse_csc)


def main():
    args = parse_args()

    print(f"Loading gene map from {args.desc}...")
    gene_map = load_gene_map(args.desc)

    print(f"Discovering gene list from {args.input}...")
    gene_list = discover_genes(args.input, args.ext, args.strand)
    if not gene_list:
        raise ValueError(f"Could not discover genes from files in {args.input}")
    print(f"  Found {len(gene_list)} genes")

    print(f"Building sparse matrix from {args.input}...")
    sparse_mtx, barcodes = build_sparse_matrix(args.input, args.ext, args.strand, gene_list)

    print(f"Writing output to {args.output}...")
    write_output(sparse_mtx, gene_list, gene_map, barcodes, args.output, args.strip_version)

    nnz = sparse_mtx.nnz
    density = nnz / (sparse_mtx.shape[0] * sparse_mtx.shape[1]) * 100 if sparse_mtx.shape[1] > 0 else 0
    print(f"Done: {sparse_mtx.shape[0]} genes x {sparse_mtx.shape[1]} cells, "
          f"{nnz} non-zero entries ({density:.2f}% density)")


if __name__ == "__main__":
    main()
