#!/usr/bin/env Rscript

suppressMessages(library(Seurat))
suppressMessages(library(stringr))
suppressMessages(library(purrr))
suppressMessages(library(Matrix))
suppressMessages(library(HDF5Array))
suppressMessages(library("optparse"))

# functions
# Define a list of options that the script accepts
option_list <- list(
  make_option(c("-o", "--Output_prefix"), type = "character", default = "custom",
              help = "Input prefix for output file, [default custom", metavar = "string"),
  make_option(c("-d", "--Input_dir"), type = "character", default = ".",
              help = "Input directory or file argument [string], [default .]", metavar = "string")
)

# Create a parser object and parse the command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

data_dir <- opt$Input_dir
oprefix <- opt$Output_prefix


# List all matrix files in the directory
mtx_files <- list.files(path = data_dir, pattern = "matrix\\.mtx$")

# Extract sample prefixes (before ".matrix.mtx")
prefixes <- sub("\\.matrix\\.mtx$", "", mtx_files)

datasets <- list()

for(i in seq_along(prefixes)) {
  prefix <- prefixes[i]
  message("Loading: ", prefix)
  
  datasets[[prefix]] <- ReadMtx(
    mtx      = paste0(prefix, ".matrix.mtx"),
    features = paste0(prefix, ".features.tsv"),
    cells    = paste0(prefix, ".barcodes.tsv")
  )
}

#-----------------------------------
# 3) Compute union of all barcodes
#-----------------------------------
all_cells <- Reduce(union, lapply(datasets, colnames))
message("Total unique cells: ", length(all_cells))

#-----------------------------------
# 4) Compute union of all genes
#-----------------------------------
all_genes <- Reduce(union, lapply(datasets, rownames))
message("Total unique genes: ", length(all_genes))

#-----------------------------------
# 5) Expand each dataset to full genes + full cells
#-----------------------------------
expanded <- lapply(datasets, function(mat) {
  # zero matrix with full genes x full cells
  out <- Matrix(0, nrow = length(all_genes), ncol = length(all_cells), sparse = TRUE)
  rownames(out) <- all_genes
  colnames(out) <- all_cells
  
  # fill in counts from the current dataset
  out[rownames(mat), colnames(mat)] <- mat
  out
})

#-----------------------------------
# 6) Combine all matrices by adding
#-----------------------------------
combined <- Reduce(`+`, expanded)

#-----------------------------------
# 7) Create final Seurat object
#-----------------------------------
merged <- CreateSeuratObject(counts = combined)

# Extract counts from the RNA assay
counts_matrix <- GetAssayData(merged, assay = "RNA", layer = "counts")

# Save as HDF5
h5_file <- paste0(oprefix, ".h5")
writeHDF5Array(counts_matrix, filepath = h5_file)