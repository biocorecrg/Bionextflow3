#!/usr/bin/env Rscript

# Concatenate bedrmod files, merging modification_names and appending data rows
# Usage: Rscript concatenate_bedrmod.R file1.bedrmod file2.bedrmod ... output.bedrmod

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: Rscript concatenate_bedrmod.R file1.bedrmod file2.bedrmod ... output.bedrmod\n")
  cat("At least 2 input files and 1 output file must be provided.\n")
  quit(status = 1)
}

# Separate input files from output file
output_file <- args[length(args)]
input_files <- args[-length(args)]

if (length(input_files) < 2) {
  cat("Error: At least 2 input files required.\n")
  quit(status = 1)
}

# Read header from first file to use as template
first_file <- input_files[1]
first_lines <- readLines(first_file)

# Find where header ends (first line not starting with #)
header_end <- which(!grepl("^#", first_lines))[1] - 1
header_lines <- first_lines[1:header_end]
column_header <- first_lines[header_end + 1]

# Extract all modification names from all files
all_modifications <- c()

for (file in input_files) {
  lines <- readLines(file)
  mod_line <- lines[grepl("^#modification_names=", lines)]

  if (length(mod_line) > 0) {
    # Extract the modification names value
    mod_names <- sub("^#modification_names=", "", mod_line)
    # Split by comma and add to all_modifications
    mods <- strsplit(mod_names, ",")[[1]]
    all_modifications <- c(all_modifications, mods)
  }
}

# Get unique modification names and sort them
unique_modifications <- unique(all_modifications)
unique_modifications <- sort(unique_modifications)
combined_mod_names <- paste(unique_modifications, collapse = ",")

# Update the modification_names line in header
header_lines <- gsub("^#modification_names=.*",
                     paste0("#modification_names=", combined_mod_names),
                     header_lines)

# Collect all data rows (non-header lines) from all files
all_data_rows <- c()

for (file in input_files) {
  lines <- readLines(file)
  # Skip all header lines (starting with #) and the column header line
  data_rows <- lines[!grepl("^#", lines) & lines != column_header & lines != ""]
  all_data_rows <- c(all_data_rows, data_rows)
}

# Write output file
output_lines <- c(header_lines, column_header, all_data_rows)
writeLines(output_lines, output_file)

cat(sprintf("Successfully concatenated %d files into %s\n",
            length(input_files), output_file))
cat(sprintf("Combined modification names: %s\n", combined_mod_names))
cat(sprintf("Total data rows: %d\n", length(all_data_rows)))
